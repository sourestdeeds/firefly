#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The backend for auto_retrieval.

@author: Steven Charles-Mindoza
"""

from ._archive import _eu, _nasa, suppress_print

from transitfit import split_lightcurve_file, run_retrieval
try:
    from lightkurve import search_lightcurve
except:
    try:
        from lightkurve import search_lightcurvefile as search_lightcurve
    except:
        raise
from traceback import format_exc
from datetime import datetime
from smtplib import SMTP_SSL
from tabulate import tabulate
from astropy.io import fits
from csv import DictWriter
from pandas import DataFrame
from shutil import rmtree, move, make_archive
import sys
import os



def _email(subject, body, to):
    username = 'transitfit.server@gmail.com'
    password = 'yvoq efzi dcib dmbm'
    sent_from = username
    # to = to ['transitfit.server@gmail.com']
    message = f'Subject: {subject}\n\n{body}'
    try:
        server = SMTP_SSL('smtp.gmail.com', 465)
        server.ehlo()
        server.login(username, password)
        server.sendmail(sent_from, to, message)
        server.close()
    except BaseException:
        # Continue on conn failure
        pass


def _TESS_filter():
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs('firefly/data', exist_ok=True)
    tess_filter_path = 'firefly/data/TESS_filter_path.csv'
    try:
        os.remove(tess_filter_path)
    except:
        pass
    tess_filter = f'{here}/data/Filters/TESS_filter.csv'
    cols = ['filter_idx', 'low_wl', 'high_wl']
    df = DataFrame(columns=cols)
    df = df.append([{'filter_idx': 0,
                     'low_wl': tess_filter}], ignore_index=True)
    df.to_csv(tess_filter_path, index=False, header=True)


def _fits(exoplanet, exo_folder, sector):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Root, Directories, Files
    files_in_dir = []
    source = f'{exo_folder}/mastDownload/TESS/'
    for r, d, f in os.walk(source):
        for item in f:
            if '.fits' in item:
                files_in_dir.append(os.path.join(r, item))
    destination = f'{exo_folder}/{exoplanet} Sector {sector}.fits'
    for files in files_in_dir:
        if files.endswith(".fits"):
            move(files, destination)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    fitsfile = destination
    with fits.open(fitsfile) as TESS_fits:
        time = TESS_fits[1].data['TIME'] + 2457000
        time += TESS_fits[1].data['TIMECORR']
        flux = TESS_fits[1].data['PDCSAP_FLUX']
        flux_err = TESS_fits[1].data['PDCSAP_FLUX_ERR']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
    write_dict = []
    for i in range(len(time)):
        write_dict.append({'Time': time[i], 'Flux': flux[i],
                           'Flux err': flux_err[i]})
    csv_name = os.path.splitext(fitsfile)[0] + '.csv'
    with open(csv_name, 'w') as f:
        columns = ['Time', 'Flux', 'Flux err']
        writer = DictWriter(f, columns)
        writer.writeheader()
        writer.writerows(write_dict)
    return fitsfile


def _MAST_query(exoplanet):
    lc = search_lightcurve(exoplanet, mission='TESS')
    if len(lc) == 0:
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    sector_list = lc .table .to_pandas()['sequence_number'] \
                     .drop_duplicates() .tolist()
    sector_list = [str(sector) for sector in sector_list]
    print(f'\nQuery from MAST returned {len(sector_list)} '
          f'data products for {exoplanet}.\n')
    lc = lc .table .to_pandas()[['observation', 
                                 'productFilename', 'size', 't_exptime']] \
            .rename(columns={'observation':'Observation'}) \
            .rename(columns={'size':'Size'}) \
            .rename(columns={'productFilename':'Product'}) 
    lc = lc[lc.t_exptime != 20].drop(['t_exptime'], axis=1)
    print(tabulate(lc, tablefmt='psql', showindex=False, headers='keys'))
    return sector_list
    

def _retrieval(exoplanet, archive='eu', curve_sample=1, nlive=300, fit_ttv=False,
               detrending_list=[['nth order', 2]],
               dynesty_sample='rslice', fitting_mode='folded',
               limb_darkening_model='quadratic', ld_fit_method='independent',
               max_batch_parameters=25, batch_overlap=2, dlogz=None, 
               maxiter=None, maxcall=None, dynesty_bounding='multi', 
               normalise=True, detrend=True):
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    exo_folder = f'firefly/{exoplanet}'
    os.makedirs(exo_folder, exist_ok=True)
    _TESS_filter()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # MAST Query
    sector_list = _MAST_query(exoplanet)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download Archive
    if archive == 'eu':
        host_T, host_z, host_r, host_logg, t0, P, t14, nan  = \
                                            _eu(exoplanet)
    elif archive == 'nasa':
        host_T, host_z, host_r, host_logg, t0, P, t14, nan = \
                                            _nasa(exoplanet)
    cols = [['t0', t0], ['P', P], ['t14', t14]]
    df = DataFrame(cols, columns=['Parameter', 'Value'])
    print('\nSplitting the lightcurve into seperate epochs'
          ' using the following parameters.\n')
    print(tabulate(df, tablefmt='psql', showindex=False, headers='keys'))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Iterate over all sectors
    curves_split, curves_delete = [], []
    for i, sector in enumerate(sector_list):
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # MAST Download
        lc = search_lightcurve(exoplanet, mission='TESS',
                                   sector=sector)
        print(f'\nDownloading MAST Lightcurve for {exoplanet} -' +
              f' TESS Sector {str(sector)}.')
        lc.download_all(download_dir=exo_folder)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Extract all light curves to a single csv file
        fitsfile = _fits(exoplanet, exo_folder, sector)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Split the Light curves
        csvfile = f'{exo_folder}/{exoplanet}_Sector_{sector}.csv'
        new_base_fname = f'sector_{sector}_split_curve'
        split_curves = split_lightcurve_file(csvfile, t0=t0, P=P, t14=t14,
                                             new_base_fname=new_base_fname)
        curves = int(curve_sample * len(split_curves))
        if curves == 0:
            curves = 1
        curves_split.append(curves)
        curves_delete.append(len(split_curves))
        print(f'\nA total of {len(split_curves)} lightcurves from TESS '
              f'Sector {sector} were created.')
        if curves == 1:
            print(f'\nA sample of {str(curves)} lightcurve from '
                  f'TESS Sector {sector} will be used.')
        else:
            print(f'\nA sample of {str(curves)} lightcurves from '
                  f'TESS Sector {sector} will be used.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths
    data_path = f'{exo_folder}/data_paths.csv'
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns=cols)
    sector_curves = dict(zip(sector_list, curves_split))
    for sector, curves in sector_curves.items():
        for i in range(curves):
            df = df.append([{'Path': f'{os.getcwd()}/{exo_folder}' +
                             f'/sector_{sector}_split_curve_{i}.csv'}],
                           ignore_index=True)
            df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
            df['Epochs'] = range(0, len(df))
    print(f'\nIn total, {len(df)} lightcurves will be fitted across all'
          ' TESS Sectors.\n')
    df.to_csv(data_path, index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = data_path
    priors = f'{exo_folder}/{exoplanet} Priors.csv'
    filters = 'firefly/data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'{exo_folder}/output_parameters'
    fitted_lightcurve_folder = f'{exo_folder}/fitted_lightcurves'
    plot_folder = f'{exo_folder}/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Run the retrieval
    run_retrieval(data, priors, filters, detrending_list=detrending_list,
                  host_T=host_T, host_logg=host_logg, host_z=host_z,
                  host_r=host_r, dynesty_sample=dynesty_sample,
                  fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                  results_output_folder=results_output_folder,
                  final_lightcurve_folder=fitted_lightcurve_folder,
                  plot_folder=plot_folder, nlive=nlive,
                  limb_darkening_model=limb_darkening_model, 
                  ld_fit_method=ld_fit_method,
                  max_batch_parameters=max_batch_parameters, 
                  batch_overlap=batch_overlap, dlogz=dlogz, 
                  maxiter=maxiter, maxcall=maxcall, 
                  dynesty_bounding=dynesty_bounding, 
                  normalise=normalise, detrend=detrend)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Cleanup
    try:
        rmtree(f'{exo_folder}/mastDownload')
        move(fitsfile, f'{exo_folder}.fits')
        os.remove(f'{exo_folder}.fits')
        os.remove(data)
        os.remove(csvfile)
        os.remove(priors)
        sector_curves = dict(zip(sector_list, curves_delete))
        for sector, curves in sector_curves.items():
            for i in range(curves):
                os.remove(f'{exo_folder}/sector_{sector}' +
                          f'_split_curve_{str(i)}.csv')
    except BaseException:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    make_archive(f'{exo_folder} {now}', format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{exoplanet}')
    rmtree(f'{exo_folder}')
    return f'{os.getcwd()}/{exo_folder} {now}.gz.tar'


def _iterable_target(exoplanet_list, archive='eu', curve_sample=1, nlive=300,
                     detrending_list=[['nth order', 2]],
                     dynesty_sample='rslice', fitting_mode='folded', fit_ttv=False,
                     limb_darkening_model='quadratic', ld_fit_method='independent',
                     max_batch_parameters=25, batch_overlap=2, dlogz=None, 
                     maxiter=None, maxcall=None, dynesty_bounding='multi', 
                     normalise=True, detrend=True, email=False, to=['transitfit.server@gmail.com'],
                     printing=False):
    for i, exoplanet in enumerate(exoplanet_list):
        try:
            # Printing suppressed within scope
            if printing == False:
                with suppress_print():
                    success = _retrieval(exoplanet, archive=archive, nlive=nlive,
                                         detrending_list=detrending_list,
                                         dynesty_sample=dynesty_sample,
                                         fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                                         limb_darkening_model=limb_darkening_model, 
                                         ld_fit_method=ld_fit_method,
                                         max_batch_parameters=max_batch_parameters, 
                                         batch_overlap=batch_overlap, dlogz=dlogz, 
                                         maxiter=maxiter, maxcall=maxcall, 
                                         dynesty_bounding=dynesty_bounding, 
                                         normalise=normalise, detrend=detrend,
                                         curve_sample=curve_sample)
            elif printing == True:
                success = _retrieval(exoplanet, archive=archive, nlive=nlive,
                                     detrending_list=detrending_list,
                                     dynesty_sample=dynesty_sample,
                                     fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                                     limb_darkening_model=limb_darkening_model, 
                                     ld_fit_method=ld_fit_method,
                                     max_batch_parameters=max_batch_parameters, 
                                     batch_overlap=batch_overlap, dlogz=dlogz, 
                                     maxiter=maxiter, maxcall=maxcall, 
                                     dynesty_bounding=dynesty_bounding, 
                                     normalise=normalise, detrend=detrend,
                                     curve_sample=curve_sample)
                print(f'\nData location: {success}\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.')
            if email == True:
                _email(f'Success: {exoplanet}',
                       f'Data location: {success} \n\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.', to=to)
        except KeyboardInterrupt:
            sys.exit('User terminated retrieval')
        except TypeError:
            trace_back = format_exc()
            if email == True:
                _email(f'Exception TypeError: {exoplanet}', trace_back, to=to)
            else:
                print(trace_back)
        except BaseException:
            trace_back = format_exc()
            if email == True:
                _email(f'Exception: {exoplanet}', trace_back, to=to)
            else:
                print(trace_back)
            pass
