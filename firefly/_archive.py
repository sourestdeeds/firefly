#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exoplanet archive tools.

@author: Steven Charles-Mindoza
"""


from transitfit import calculate_logg
from datetime import datetime, timedelta
from tabulate import tabulate
from pandas import DataFrame, read_csv
from fuzzywuzzy import process
from natsort import natsorted
import numpy as np
import random
import sys
import os



class suppress_print():
    def __enter__(self):
        self.original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self.original_stdout


def _load_csv():
    _download_nasa()
    here = os.path.dirname(os.path.abspath(__file__))
    nasa_csv = 'firefly/data/nasa.csv.gz'
    mast_csv = f'{here}/data/Filters/MAST_lc.csv.xz'
    global exo, mast
    exo = read_csv(nasa_csv)
    mast = read_csv(mast_csv)        
    

def _search(exoplanet):
    exo_list = exo[['pl_name', 'tic_id']] \
              .dropna() .drop_duplicates('pl_name') \
              .drop(['tic_id'], axis=1) .values .tolist()
    exo_list = [j for i in exo_list for j in i]
    exo_list = natsorted(exo_list)
    
    ratios = process.extract(exoplanet,exo_list)
    highest = process.extractOne(exoplanet,exo_list)
    return highest, ratios


def _tic(exoplanet):
    '''
    Returns a dataframe of all planet names and tic ID's.

    Parameters
    ----------
    None

    Returns
    -------
    TIC ID of an exoplanet target.

    '''
    df = exo[['pl_name','tic_id','soltype','rowupdate','disc_facility']] \
        .sort_values('pl_name') \
        .drop_duplicates('pl_name', keep='last')
    tic_df = df .drop(['soltype', 'rowupdate', 'disc_facility'], axis=1).dropna() \
                .set_index('pl_name')
    target = tic_df.loc[[exoplanet]]
    tic_id = target['tic_id'] .values .tolist()[0]
    return tic_id


def _lc(exoplanet):
    '''
    # https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html
    out = DataFrame()
    for i in range(1,32):
        a = read_csv(f'TESS_curl/tesscurl_sector_{str(i)}_lc.sh')
        a['links'] = a['#!/bin/sh'].str.split().str[-1]
        a = a.drop(['#!/bin/sh'], axis=1)
        out = out.append(a)
    out.to_csv('MAST_lc.csv.gz', index=False)
    '''
    tic_id = _tic(exoplanet).replace('TIC ', '')
    lc_links = mast[mast['links'].str.contains(tic_id)] .values .tolist()
    lc_list = [i for j in lc_links for i in j]
    lc_test = [int(i[-30:-15]) for j in lc_links for i in j]
    lc_links = []
    for i in range(len(lc_test)):
        if int(tic_id) == lc_test[i]:
            lc_links.append(lc_list[i])
    return lc_links, tic_id
    

def _download_nasa():
    os.makedirs('firefly/data', exist_ok=True)
    download_link =  \
        'https://exoplanetarchive.ipac.caltech.edu/' +\
        'TAP/sync?query=select+' +\
        'pl_name,tic_id,pl_orbper,pl_orbsmax,pl_radj,pl_orbeccen,ttv_flag,' +\
        'st_teff,st_rad,st_mass,st_met,st_logg,pl_tranmid,pl_trandur,' +\
        'pl_orbincl,pl_orblper,soltype,rowupdate,disc_facility' +\
        '+from+ps&format=csv'
    nasa_csv = 'firefly/data/nasa.csv.gz'
    if not os.path.exists(nasa_csv):
        print('Caching the NASA Exoplanet Archive.')
        df = read_csv(download_link)
        df.to_csv(nasa_csv, index=False)
    ten_days_ago = datetime.now() - timedelta(days=10)
    filetime = datetime.fromtimestamp(os.path.getctime(nasa_csv))
    if filetime < ten_days_ago:
        print('NASA Archive is 10 days old. Updating.')
        df = read_csv(download_link)
        df.to_csv(nasa_csv, index=False)
    else:
        pass
    return nasa_csv


def _nasa_full():
    '''
    Downloads the entire NASA Exoplanet Archive

    Returns
    -------
    data/nasa_full.csv

    '''
    download_link =  \
    'https://exoplanetarchive.ipac.caltech.edu/' +\
    'TAP/sync?query=select+*+from+ps&format=csv'
    df = read_csv(download_link)
    df.to_csv('firefly/data/nasa_full.csv.xz', index=False)


def estimate_t14(Rp, Rs, a, P):
    '''
    Estimates t14 in minutes, if P is in days
    '''
    AU = 1.495978707e11
    R_sun = 6.957e8
    R_jup = 71492000
    return (Rp * R_jup + Rs * R_sun)/(np.pi * a * AU) * P * 24 * 60


def _IQR(df):
    # Only keep IQR of data
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    trueList = ~((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR)))
    return df[trueList]
  
    
def _nasa(exoplanet, save=True):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download NASA archive
    _download_nasa()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Read in nasa.csv
    exo_archive = exo.set_index('pl_name')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    try:
        df = DataFrame(exo_archive).loc[[exoplanet]]
    except KeyError:
        sys.exit('The chosen target is either spelt incorrectly, or does not '
                 'exist in the NASA archive.')
    tic = df['tic_id'] .drop_duplicates() .max()
    # Fix t0 on most recent centred transit
    t0 = df['pl_tranmid'] .max()
    # Only keep IQR of data
    df = _IQR(df)
    # Average the rest
    s = df.mean()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Values for calculation
    P = s.loc['pl_orbper']
    a = s.loc['pl_orbsmax']
    i = s.loc['pl_orbincl']
    w = s.loc['pl_orblper']
    ecc = s.loc['pl_orbeccen']
    rp = s.loc['pl_radj']
    rs = s.loc['st_rad']
    z = s.loc['st_met']
    ms = s.loc['st_mass']
    logg = s.loc['st_logg']
    T = s.loc['st_teff']
    t14 = s.loc['pl_trandur'] * 60
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    host_T = (T, T * 1e-2)
    host_z = (z, np.abs(z * 1e-2))
    host_r = (rs, rs * 1e-2)
    host_logg = (logg, logg * 1e-2)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Nan checks
    # G = 6.67408e-11
    # AU = 1.495978707e11
    # sol = 1.98847e30
    # P
    # if (np.isnan(P) and not np.isnan(a)):
    #     # Added small correction factor to star mass
    #     P = 2 * np.pi * np.sqrt((a * AU)**3 / 
    #                             (G * 0.963 * ms * sol)) / (60 * 60 * 24)
    # # a
    # elif (np.isnan(a) and not np.isnan(P)):
    #     a = (((P * 24 * 60 * 60)**2 * G * ms * sol / 
    #           (4 * np.pi**2))**(1 / 3)) / AU
    # w
    if np.isnan(w):
        w = 90
    # ecc
    elif np.isnan(ecc):
        ecc = 0
    # t14
    elif np.isnan(t14):
        t14 = estimate_t14(rp, rs, a, P)
    # logg
    elif np.isnan(logg):
        logg, err_logg = calculate_logg((ms, ms * 1e-2), (rs, rs * 1e-2))
        host_logg = (logg, err_logg)
    # z
    elif (np.isnan(z) or z==0):
        host_z = (0, 0.1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Exoplanet Priors to TransitFit
    radius_const = 0.1027626851
    cols = [['P', 'gaussian', P, P * 1e-4, ''],
            ['t0', 'gaussian', t0, 7e-3, ''],
            ['a', 'gaussian', a, a * 1e-2, ''],
            ['inc', 'gaussian', i, i * 1e-2, ''],
            ['w', ['fixed' if (w==90 or w==0) else 'gaussian'][0], w,
             ['' if (w==90 or w==0) else np.abs(w * 1e-2)][0], ''],
            ['ecc', ['fixed' if ecc==0 else 'gaussian'][0], ecc,
             ['' if ecc==0 else ecc * 1e-2][0], ''],
            ['rp', 'uniform',
             0.9 * radius_const * rp / rs,
             1.1 * radius_const * rp / rs, 0]]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save the priors
    priors_csv = DataFrame(cols, columns=['Parameter', 'Distribution',
                                          'Input_A', 'Input_B', 'Filter'])
    if save == True:
        priors = f'firefly/{exoplanet}/{exoplanet} Priors.csv'
        if not os.path.exists(priors):
            priors_csv.to_csv(priors, index=False, header=True)
        else:
            pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # For printing variables only
    cols = [['P', 'gaussian', P, P * 1e-4, ''],
            ['t0', 'gaussian', t0, 7e-3, ''],
            ['a', 'gaussian', a, a * 1e-2, ''],
            ['inc', 'gaussian', i, i * 1e-2, ''],
            ['w', ['fixed' if (w==90 or w==0) else 'gaussian'][0], w,
             ['' if (w==90 or w==0) else np.abs(w * 1e-2)][0], ''],
            ['ecc', ['fixed' if ecc==0 else 'gaussian'][0], ecc,
             ['' if ecc==0 else ecc * 1e-2][0], ''],
            ['rp', 'uniform',
             0.9 * radius_const * rp / rs,
             1.1 * radius_const * rp / rs, 0],
            ['host_T', 'fixed', host_T[0], host_T[1], ''],
            ['host_z', 'fixed', host_z[0], host_z[1], ''],
            ['host_r', 'fixed', host_r[0], host_r[1], ''],
            ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input A', 'Input B', 'Filter'])
    nan = repack.isnull().values.any()
    print(f'\nPriors generated from the NASA Archive for {exoplanet}'
          f' ({tic}).\n')
    print(tabulate(repack, tablefmt='psql', showindex=False, headers='keys'))
    return host_T, host_z, host_r, host_logg, t0, P, t14, nan, repack


def tess_viable(k=10, survey=None):
    '''
    Currently there are 377 tess targets with full prior and lightcurve sets.

    Parameters
    ----------
    k : int, optional
        Returns a random set of targets. The default is 10.
    survey : str, optional
        Choice of survey to filter by. The default is None.

    Returns
    -------
    targets : str
        Viable tess target list.

    '''
    here = os.path.dirname(os.path.abspath(__file__))
    tess = f'{here}/data/Filters/tess_viable.csv'
    tess_ttv = f'{here}/data/Filters/tess_ttv_viable.csv'
    targets = read_csv(tess)['Exoplanet'] 
    ttv_targets = read_csv(tess_ttv)['Exoplanet']
    if survey != None:
        targets = [s for s in targets if survey in s]
        ttv_targets = [s for s in ttv_targets if survey in s]
    all_targets = natsorted(targets)
    ttv_targets = natsorted(ttv_targets)
    try:
        targets = random.sample(targets, k=k)
    except ValueError:
        k = len(targets)
        targets = random.sample(targets, k=k)
    return targets, all_targets, ttv_targets


def generate_tess_viable():
    _download_nasa()
    here = os.path.dirname(os.path.abspath(__file__))
    nasa_csv = 'firefly/data/nasa.csv.gz'
    mast_csv = f'{here}/data/Filters/MAST_lc.csv.xz'
    global exo, mast
    exo = read_csv(nasa_csv)
    mast = read_csv(mast_csv)
    
    exo_list = exo[['pl_name', 'tic_id']] \
              .dropna() .drop_duplicates('pl_name') \
              .drop(['tic_id'], axis=1) .values .tolist()
    exo_list = [j for i in exo_list for j in i]
    
    viable = []
    products = []
    for i, exoplanet in enumerate(exo_list):
        lc_links, tic_id = _lc(exoplanet)
        nan = _check_nan(exoplanet)
        if (len(lc_links)==0 or nan==True):
            pass
        else:
            print(f'{exoplanet} viable.')
            viable.append(exoplanet)
            products.append(len(lc_links))
    data = {'Exoplanet':viable, 'Products':products}
    df = DataFrame(data).sort_values('Exoplanet')
    df.to_csv(f'{here}/data/Filters/tess_viable.csv', index=False)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # TTV
    ttv_list = exo[['pl_name', 'tic_id', 'ttv_flag']] 
    ttv_list = ttv_list[ttv_list!=0] . dropna() \
              .drop_duplicates('pl_name') \
              .drop(['tic_id', 'ttv_flag'], axis=1) .values .tolist()
    ttv_list = [j for i in ttv_list for j in i]
    ttv_list = natsorted(ttv_list)
    # ttv_list = [s for s in ttv_list if 'Kepler' not in s]
    # ttv_list = [s for s in ttv_list if 'K2' not in s]
    # ttv_list = [s for s in ttv_list if 'KOI' not in s]
    viable_ttv = []
    products_ttv = []
    for i, exoplanet in enumerate(ttv_list):
        lc_links, tic_id = _lc(exoplanet)
        nan = _check_nan(exoplanet)
        if (len(lc_links)==0 or nan==True):
            pass
        else:
            print(f'{exoplanet} viable.')
            viable_ttv.append(exoplanet)
            products_ttv.append(len(lc_links))            
    data = {'Exoplanet':viable_ttv, 'Products':products_ttv}
    df = DataFrame(data).sort_values('Exoplanet')
    df.to_csv(f'{here}/data/Filters/tess_ttv_viable.csv', index=False)
    



def _check_nan(exoplanet, printing=False):
    if printing == False:
        with suppress_print():
            nan = _nasa(exoplanet, save=False)
            nan = nan[7]
    elif printing == True:
        nan = _nasa(exoplanet, save=False)
        nan = nan[7]
    return nan
