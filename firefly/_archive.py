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
import numpy as np
import sys
import os



class suppress_print():
    def __enter__(self):
        self.original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self.original_stdout     
        
    

def _download_nasa():
    os.makedirs('firefly/data', exist_ok=True)
    download_link =  \
        'https://exoplanetarchive.ipac.caltech.edu/' +\
        'TAP/sync?query=select+' +\
        'pl_name,tic_id,pl_orbper,pl_orbsmax,' +\
        'pl_radj,pl_orbeccen,' +\
        'st_teff,st_tefferr1,st_rad,st_raderr1,st_mass,st_masserr1,' +\
        'st_met,st_meterr1,st_logg,st_loggerr1,pl_tranmid,pl_trandur,' +\
        'pl_orbincl,pl_orblper,soltype,rowupdate,disc_facility' +\
        '+from+ps&format=csv'
    nasa_csv = 'firefly/data/nasa.csv.gz'
    if not os.path.exists(nasa_csv):
        print('nasa.csv does not exist, downloading.')
        df = read_csv(download_link)
        df.to_csv(nasa_csv, index=False)
    ten_days_ago = datetime.now() - timedelta(days=10)
    filetime = datetime.fromtimestamp(os.path.getctime(nasa_csv))
    if filetime < ten_days_ago:
        print('nasa.csv is 10 days old, updating.')
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
    nasa_csv = _download_nasa()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Read in nasa.csv
    exo_archive = read_csv(nasa_csv, index_col='pl_name') 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    try:
        df = DataFrame(exo_archive).loc[[exoplanet]]
    except KeyError:
        sys.exit('The chosen target is either spelt incorrectly, or does not '
                 'exist in the NASA archive.')
    tic = df['tic_id'] .drop_duplicates() .values .tolist()[0]
    # Only keep IQR of data
    df = _IQR(df)
    # Fix t0 on first centred transit
    try:
        t0 = df['pl_tranmid'] .dropna()[0]
    except IndexError:
        sys.exit('No value for t0 was found in the NASA Archive'
                 f' for {exoplanet}.')
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
    m_s = s.loc['st_mass']  
    logg = s.loc['st_logg']
    T = s.loc['st_teff']
    t14 = s.loc['pl_trandur'] * 60
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    host_T = (s.loc['st_teff'], s.loc['st_tefferr1'])
    host_z = (s.loc['st_met'], s.loc['st_meterr1'])
    host_r = (s.loc['st_rad'], s.loc['st_raderr1'])
    host_logg = (s.loc['st_logg'], s.loc['st_loggerr1'])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Nan checks
    if np.isnan(t14):
        t14 = estimate_t14(rp, rs, a, P)
    elif np.isnan(logg):
        logg, err_logg = calculate_logg((s.loc['st_mass'],
                                     s.loc['st_masserr1']),
                                    (s.loc['st_rad'],
                                     s.loc['st_raderr1']))
        host_logg = (logg, err_logg)
    elif np.isnan(z):
        host_z = (0, 0.1)
    elif np.isnan(w):
        w = 90
    elif np.isnan(ecc):
        ecc = 0
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Exoplanet Priors to TransitFit
    radius_const = 0.1027626851
    cols = [['P', 'gaussian', P, P * 1e-5 ''],
            ['t0', 'gaussian', t0, 7e-3, ''],
            ['a', 'gaussian', a, a * 0.1, ''],
            ['inc', 'gaussian', i, i * 0.1, ''],
            ['w', ['fixed' if w==90 else 'gaussian'][0], w, 
             ['' if w==90 else w * 0.1][0], ''],
            ['ecc', ['fixed' if ecc==0 else 'gaussian'][0], ecc, 
             ['' if ecc==0 else ecc * 0.1][0], ''],
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
    cols = [['P', 'gaussian', P, P * 1e-5, ''],
            ['t0', 'gaussian', t0, 7e-3, ''],
            ['a', 'gaussian', a, a * 0.1, ''],
            ['inc', 'gaussian', i, i * 0.1, ''],
            ['w', ['fixed' if w==90 else 'gaussian'][0], w, 
             ['' if w==90 else w * 0.1][0], ''],
            ['ecc', ['fixed' if ecc==0 else 'gaussian'][0], ecc, 
             ['' if ecc==0 else ecc * 0.1][0], ''],
            ['rp', 'uniform',
             0.9 * radius_const * rp / rs,
             1.1 * radius_const * rp / rs, 0],
            ['host_T', 'fixed', host_T[0], host_T[1], ''],
            ['host_z', 'fixed', host_z[0], host_z[1], ''],
            ['host_r', 'fixed', host_r[0], host_r[1], ''],
            ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]           
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input_A', 'Input_B', 'Filter'])
    nan = repack.isnull().values.any()
    print(f'\nPriors generated from the NASA Archive for {exoplanet}'
          f' ({tic}).\n')
    print(tabulate(repack, tablefmt='psql', showindex=False, headers='keys'))
    return host_T, host_z, host_r, host_logg, t0, P, t14, nan, repack


def _check_nan(exoplanet, printing=False):
    if printing == False:
        with suppress_print():
            nan = _nasa(exoplanet, save=False)
            nan = nan[7]
    elif printing == True:
        nan = _nasa(exoplanet, save=False)
        nan = nan[7]
    return nan
