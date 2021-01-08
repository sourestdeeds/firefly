#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The backend for auto_retrieval.

@author: Steven Charles-Mindoza
"""

from ._archive import _eu, _nasa
from .query import _lc
from transitfit import split_lightcurve_file, run_retrieval
from datetime import datetime
from smtplib import SMTP_SSL
from tabulate import tabulate
from astropy.io import fits
from csv import DictWriter
from pandas import DataFrame, read_csv
from shutil import rmtree, make_archive
from email.message import EmailMessage
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
from math import ceil
import sys
import os
import logging
import random

logging.getLogger('googleapiclient.discovery_cache').setLevel(logging.ERROR)

def _gdrive(archive_name, fitting_mode):
    gauth = GoogleAuth()
    here = os.path.dirname(os.path.abspath(__file__))
    gauth.LoadCredentialsFile(f'{here}/data/Filters/fireflycreds.txt')
    # gauth.LoadCredentialsFile('fireflycreds.txt')
    # gauth.LocalWebserverAuth()
    # gauth.SaveCredentialsFile('fireflycreds.txt')
    drive = GoogleDrive(gauth)
    if fitting_mode == 'all':
        folder = drive \
                 .ListFile({'q': "title='all' and trashed=False"}) \
                 .GetList()[0]
    elif fitting_mode == 'folded':
        folder = drive \
                 .ListFile({'q': "title='folded' and trashed=False"}) \
                 .GetList()[0]
    elif fitting_mode == 'batched':
        folder = drive \
                 .ListFile({'q': "title='batched' and trashed=False"}) \
                 .GetList()[0]
    elif fitting_mode == '2_stage':
        folder = drive \
                 .ListFile({'q': "title='2_stage' and trashed=False"}) \
                 .GetList()[0]
    print('\nUploading results to googledrive.')
    file = drive.CreateFile({'title':f'{archive_name}.tar.gz', 
                             'parents':[{'id':folder['id']}]})
    file.SetContentFile(f'firefly/{archive_name}.tar.gz')
    file.Upload()
    
    keys = file.keys()
    if file['shared']:
        link = 'https://drive.google.com/drive/folder/' + file['id'] +\
               '?usp=sharing'
    elif 'webContentLink' in keys:
        link = file['webContentLink']
    elif 'webViewLink' in keys:
        link = file['webViewLink']
    else:
        link = 'No link available.'
    
    if 'name' in keys:
        name = file['name']
    else:
        name = file['id']
    return link


def _email(subject, body, to):
    username = 'transitfit.server@gmail.com'
    password = 'spearnet1!'
    message = f'Subject: {subject}\n\n{body}'
    
    email = EmailMessage()
    email['Subject'] = subject
    email['From'] = 'transitfit.server@gmail.com'
    email['To'] = to
    email.set_content(body, subtype='html')
    #email.add_alternative(args, kw)
    #email.add_attachment()
    
    try:
        with SMTP_SSL('smtp.gmail.com', 465) as s:
            #s.ehlo()
            s.login(username, password)
            s.send_message(email)
        
        #server.sendmail(sent_from, to, message)
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


def _fits(exoplanet, exo_folder, cache):
    print(f'\nSearching MAST for {exoplanet}.')
    lc_links, tic_id = _lc(exoplanet)
    if len(lc_links) == 0:
        rmtree(exo_folder)
        print(f'Search result contains no data products for {exoplanet}.')
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    print(f'\nQuery from MAST returned {len(lc_links)} '
          f'data products for {exoplanet} (TIC {tic_id}).\n')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    csv_in_dir = []
    fitsname = []
    sector_list = []
    for j, fitsfile in enumerate(lc_links):
        with fits.open(fitsfile, cache=cache) as TESS_fits:
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
        mast_name = fitsfile[-55:].replace('.fits', '')
        os.makedirs(f'{source}/{mast_name}', exist_ok=True)
        csv_name = f'{source}/{mast_name}/{mast_name}.csv'
        with open(csv_name, 'w') as f:
            columns = ['Time', 'Flux', 'Flux err']
            writer = DictWriter(f, columns)
            writer.writeheader()
            writer.writerows(write_dict)
        csv_in_dir.append(f'{os.getcwd()}/{csv_name}')
        fitsname.append(fitsfile[-55:])
        sector_list.append(int(fitsfile[-36:-32]))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Print search result 
    tess_sector = [f'TESS Sector {str(sector)}' for sector in sector_list]
    _ = {'Sector':tess_sector, 'Product':fitsname}
    df = DataFrame(_)
    print()
    print(tabulate(df, tablefmt='psql', showindex=False, headers='keys'))
    print('\nSplitting up the lightcurves into seperate epochs:\n')
    # csv_in_dir = []
    # for r, d, f in os.walk(source):
    #     for item in f:
    #         if '.csv' in item:
    #             csv_in_dir.append(os.path.join(r, item))
    return csv_in_dir
    

def _retrieval(
        # Firefly Interface
        exoplanet, 
        archive='eu', 
        curve_sample=1, 
        clean=False,
        cache=False,
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
        host_T, host_z, host_r, host_logg, t0, P, t14, nan, repack  = \
                                            _eu(exoplanet)
    elif archive == 'nasa':
        try:
            host_T, host_z, host_r, host_logg, t0, P, t14, nan, repack = \
                                                _nasa(exoplanet)
        except Exception:
            os.remove('firefly/data/nasa.csv.gz')
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
    csv_in_dir = _fits(exoplanet, exo_folder, cache)
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
    results = read_csv(f'{exo_folder}/output_parameters/Complete_results.csv')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Cleanup
    if clean == True:
        rmtree(f'{exo_folder}/mastDownload')
        os.remove(data)
        os.remove(priors)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    archive_name = f'{exoplanet} {now}'
    archive_folder = f'firefly/{archive_name}'
    make_archive(archive_folder, format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{exoplanet}')
    rmtree(f'{exo_folder}')
    return archive_name, repack, results

