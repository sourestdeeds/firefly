#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The backend for auto_retrieval.

@author: Steven Charles-Mindoza
"""

from ._archive import _eu, _nasa
from .query import tic, _lc

from transitfit import split_lightcurve_file, run_retrieval
try:
    from lightkurve import search_lightcurve
except:
    try:
        from lightkurve import search_lightcurvefile as search_lightcurve
    except:
        raise
from datetime import datetime
from smtplib import SMTP_SSL
from tabulate import tabulate
from astropy.io import fits
from csv import DictWriter
from pandas import DataFrame
from shutil import rmtree, make_archive
from math import ceil
import sys
import os
import random



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


def _fits(exoplanet, exo_folder, clean):
    tic_id = tic(exoplanet).replace('TIC ', '')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Path to downloaded fits
    fits_in_dir = []
    source = f'{exo_folder}/mastDownload/TESS/'
    for r, d, f in os.walk(source):
        for item in f:
            # Delete the fast-lc files
            if ('fast-lc.fits') in item:
                print('Deleting:', *f)
                rmtree(r)
            # Check we have the correct TIC
            elif not tic_id in item:
                print('Deleting:', *f)
                rmtree(r)
            # Keep the _lc files
            elif tic_id in item:
                fits_in_dir.append(os.path.join(r, item))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    for i, fitsfile in enumerate(fits_in_dir):
        with fits.open(fitsfile) as TESS_fits:
            time = TESS_fits[1].data['TIME'] + 2457000
            time += TESS_fits[1].data['TIMECORR']
            flux = TESS_fits[1].data['PDCSAP_FLUX']
            flux_err = TESS_fits[1].data['PDCSAP_FLUX_ERR']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
        if clean == True:
            os.remove(fitsfile)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Path to downloaded csv   
    csv_in_dir = []
    for r, d, f in os.walk(source):
        for item in f:
            if '.csv' in item:
                csv_in_dir.append(os.path.join(r, item))
    return csv_in_dir



def _fits_quick(exoplanet, exo_folder, clean):
    lc_links, tic_id = _lc(exoplanet)
    print(f'\nSearching MAST for {exoplanet} ({tic_id}).')
    if len(lc_links) == 0:
        rmtree(exo_folder)
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    print(f'\nQuery from MAST returned {len(lc_links)} '
          f'data products for {exoplanet} ({tic_id}).')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    for j, fitsfile in enumerate(lc_links):
        with fits.open(fitsfile) as TESS_fits:
            time = TESS_fits[1].data['TIME'] + 2457000
            time += TESS_fits[1].data['TIMECORR']
            flux = TESS_fits[1].data['PDCSAP_FLUX']
            flux_err = TESS_fits[1].data['PDCSAP_FLUX_ERR']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
        write_dict = []
        for i in range(len(time)):
            write_dict.append({'Time': time[i], 'Flux': flux[i],
                               'Flux err': flux_err[i]})
        source = f'{exo_folder}/mastDownload'
        os.makedirs(source+f'/{str(j)}', exist_ok=True)
        csv_name = f'{exo_folder}/mastDownload/{str(j)}/reduced_{str(j)}.csv'
        with open(csv_name, 'w') as f:
            columns = ['Time', 'Flux', 'Flux err']
            writer = DictWriter(f, columns)
            writer.writeheader()
            writer.writerows(write_dict)
        if clean == True:
            os.remove(fitsfile)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Path to downloaded csv   
    csv_in_dir = []
    for r, d, f in os.walk(source):
        for item in f:
            if '.csv' in item:
                csv_in_dir.append(os.path.join(r, item))
    return csv_in_dir


def _MAST_query(exoplanet, exo_folder):
    tic_id = tic(exoplanet)
    print(f'\nSearching MAST for {exoplanet} ({tic_id}).')
    lc = search_lightcurve(exoplanet, mission='TESS', radius=0.01)
    if len(lc) == 0:
        rmtree(exo_folder)
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    sector_list = lc .table .to_pandas()['sequence_number'] \
                      .drop_duplicates() .tolist()
    sector_list = [str(sector) for sector in sector_list]
    lc = lc .table .to_pandas()[['target_name', 'observation', 
                                 'productFilename', 't_exptime']] \
            .rename(columns={'target_name':'Target'}) \
            .rename(columns={'observation':'Observation'}) \
            .rename(columns={'productFilename':'Product'}) 
    lc = lc[lc.t_exptime != 20].drop(['t_exptime'], axis=1)
    print(f'\nQuery from MAST returned {len(lc)} '
          f'data products for {exoplanet} ({tic_id}).\n')
    print(tabulate(lc, tablefmt='psql', showindex=False, headers='keys'))
    return sector_list
    

def _retrieval(
        # Firefly Interface
        exoplanet, 
        archive='eu', 
        curve_sample=1, 
        clean=False,
        # TransitFit Variables
        cutoff=0.25, 
        window=2.5,
        nlive=300, 
        fit_ttv=False, 
        detrending_list=[['nth order', 2]],
        dynesty_sample='auto', 
        fitting_mode='auto',
        limb_darkening_model='quadratic', 
        ld_fit_method='independent',
        max_batch_parameters=25, 
        batch_overlap=2, 
        dlogz=None, 
        maxiter=None, 
        maxcall=None, 
        dynesty_bounding='multi', 
        normalise=True, 
        detrend=True
):
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    exo_folder = f'firefly/{exoplanet}'
    os.makedirs(exo_folder, exist_ok=True)
    _TESS_filter()
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
    # Split the Light curves
    split_curve_in_dir = []
    csv_in_dir = _fits_quick(exoplanet, exo_folder, clean)
    for i, csvfile in enumerate(csv_in_dir):
        split_curves = split_lightcurve_file(csvfile, t0=t0, P=P, t14=t14, 
                                             cutoff=cutoff, window=window)
        split_curves = [s + '.csv' for s in split_curves]
        split_curve_in_dir.append(split_curves)
        if clean == True:
            os.remove(csvfile)
    split_curve_in_dir = [i for sub in split_curve_in_dir for i in sub]
    print(f'\nA total of {len(split_curve_in_dir)} lightcurves '
          'were generated.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Sort the files into ascending order and take random sample
    curves = ceil(curve_sample * len(split_curve_in_dir))
    split_curve_in_dir = random.sample(split_curve_in_dir, k=int(curves))
    split_curve_in_dir.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths     
    data_path = f'{exo_folder}/data_paths.csv'
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns=cols)
    for i, split_curve in enumerate(split_curve_in_dir):
        df = df.append([{'Path': split_curve}], ignore_index=True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    print(f'\nA random sample of {len(df)} lightcurves will be fitted'
          ' across all TESS Sectors.\n')
    df.to_csv(data_path, index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = data_path
    priors = f'{exo_folder}/{exoplanet} Priors.csv'
    filters = 'firefly/data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'{exo_folder}/output_parameters'
    fitted_lightcurve_folder = f'{exo_folder}/fitted_lightcurves'
    plot_folder = f'{exo_folder}/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Run the retrieval
    run_retrieval(
            data, 
            priors, 
            filters, 
            host_T=host_T, 
            host_logg=host_logg, 
            host_z=host_z, 
            host_r=host_r,
            nlive=nlive,
            fit_ttv=fit_ttv,
            detrending_list=detrending_list,
            dynesty_sample=dynesty_sample,
            fitting_mode=fitting_mode, 
            limb_darkening_model=limb_darkening_model, 
            ld_fit_method=ld_fit_method,
            max_batch_parameters=max_batch_parameters, 
            batch_overlap=batch_overlap, 
            dlogz=dlogz, 
            maxiter=maxiter, 
            maxcall=maxcall, 
            dynesty_bounding=dynesty_bounding, 
            normalise=normalise, 
            detrend=detrend,
            results_output_folder=results_output_folder,
            final_lightcurve_folder=fitted_lightcurve_folder,
            plot_folder=plot_folder
        )
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Cleanup
    if clean == True:
        rmtree(f'{exo_folder}/mastDownload')
        os.remove(data)
        os.remove(priors)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    make_archive(f'{exo_folder} {now}', format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{exoplanet}')
    rmtree(f'{exo_folder}')
    return f'{os.getcwd()}/{exo_folder} {now}.gz.tar'

