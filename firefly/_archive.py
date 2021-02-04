#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exoplanet archive tools.

@author: Steven Charles-Mindoza
"""


from transitfit import calculate_logg
from datetime import datetime, timedelta
from tabulate import tabulate
from pandas import DataFrame, read_csv, Categorical, concat
from fuzzywuzzy import process
from natsort import natsorted
from math import ceil
import numpy as np
import random
import sys
import os


class NaNError(Exception):
    pass

class suppress_print():
    def __enter__(self):
        self.original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self.original_stdout


def _load_csv():
    '''
    Coverts all 4 archives to have the same column structure.

    '''
    _download_archive()
    here = os.path.dirname(os.path.abspath(__file__))
    nasa_csv = 'firefly/data/nasa.csv.gz'
    eu_csv = 'firefly/data/eu.csv.gz'
    oec_csv = 'firefly/data/oec.csv.gz'
    org_csv = 'firefly/data/org.csv.gz'
    mast_csv = f'{here}/data/Filters/MAST_lc.csv.xz'
    global exo_nasa, exo_eu, exo_oec, exo_org, mast
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # NASA
    exo_nasa = read_csv(nasa_csv)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # EU
    col_subset_eu = ['# name', 'orbital_period', 'semi_major_axis', 'radius',
                     'eccentricity', 'inclination', 'tzero_tr', 'star_teff',
                     'star_teff_error_max', 'star_radius',
                     'star_radius_error_max', 'star_mass',
                     'star_metallicity', 'star_metallicity_error_max', 'omega']
    exo_eu = read_csv(eu_csv, usecols=col_subset_eu)
    exo_eu.columns = ['pl_name', 'pl_radj', 'pl_orbper', 'pl_orbsmax',
                      'pl_orbeccen', 'pl_orbincl', 'pl_orblper', 'pl_tranmid', 'st_met',
                      'st_meterr1', 'st_mass','st_rad',
                      'st_raderr1', 'st_teff', 'st_tefferr1']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # OEC
    col_subset_oec = ['name', 'radius', 'period',
                     'semimajoraxis', 'eccentricity', 'periastron', 'inclination',
                     'hoststar_mass', 'hoststar_radius',
                     'hoststar_metallicity', 'hoststar_temperature']
    exo_oec = read_csv(oec_csv, usecols=col_subset_oec)
    exo_oec.columns = ['pl_name', 'pl_radj', 'pl_orbper', 'pl_orbsmax',
                      'pl_orbeccen', 'pl_orblper','pl_orbincl', 'st_mass',
                      'st_rad', 'st_met', 'st_teff']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # ORG
    col_subset_org = ['NAME', 'R', 'PER',
                     'SEP', 'ECC', 'I', 'TT',
                     'MSTAR', 'RSTAR', 'RSTARUPPER',
                     'FE', 'FEUPPER', 'LOGG', 'LOGGUPPER',
                     'TEFF', 'TEFFUPPER']
    exo_org = read_csv(org_csv, usecols=col_subset_org)
    exo_org.columns = ['pl_orbeccen', 'st_met', 'st_meterr1', 'pl_orbincl',
                      'st_logg', 'st_loggerr1','st_mass', 'pl_name',
                      'pl_orbper', 'pl_radj', 'st_rad', 'st_raderr1', 'pl_orbsmax',
                      'st_teff', 'st_tefferr1', 'pl_tranmid']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # MAST
    mast = read_csv(mast_csv)
    

def _search(exoplanet):
    exo_list = exo_nasa[['pl_name', 'tic_id']] \
              .dropna() \
              .drop(['tic_id'], axis=1) .values .tolist()
    exo_list = [j for i in exo_list for j in i]
    exo_list = natsorted(exo_list)
    
    ratios = process.extract(exoplanet,exo_list)
    highest = process.extractOne(exoplanet,exo_list)
    return highest, ratios


def _search_all(exoplanet):
    archive_list = [exo_nasa, exo_eu, exo_oec, exo_org]
    exo_archive = concat(archive_list)# .set_index('pl_name')
    exo_list = exo_archive[['pl_name']] \
              .dropna() .drop_duplicates('pl_name') \
              .values .tolist()
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
    df = exo_nasa[['pl_name','tic_id']] \
        .sort_values('pl_name') \
        .drop_duplicates('pl_name', keep='last')
    tic_df = df .dropna() .set_index('pl_name')
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


def tess_candidate(count=18):
    '''
    Finds all candidate links with 18 sector observations.

    '''
    df = mast.sort_values('links').reset_index()
    df = df.groupby(df.links.str[-30:-15].astype('int32')).size().reset_index(name='Products')
    df = df.loc[df['Products']==count]
    df = df.rename(columns={'links':'TIC ID'})
    df.to_csv('firefly/data/{count}_sector_candidates.csv', index=False)
    

def _download_archive():
    os.makedirs('firefly/data', exist_ok=True)
    download_links = [ \
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # NASA
        # 'https://exoplanetarchive.ipac.caltech.edu/' +\
        # 'TAP/sync?query=select+' +\
        # 'pl_name,tic_id,pl_orbper,pl_orbsmax,pl_radj,pl_orbeccen,ttv_flag,' +\
        # 'st_teff,st_rad,st_mass,st_met,st_logg,pl_tranmid,pl_trandur,' +\
        # 'st_tefferr1,st_raderr1,st_meterr1,st_loggerr1,' +\
        # 'pl_orbincl,pl_orblper' +\
        # '+from+ps&format=csv',
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # NASA COMP
        'https://exoplanetarchive.ipac.caltech.edu/' +\
        'TAP/sync?query=select+' +\
        'pl_name,tic_id,pl_orbper,pl_orbsmax,pl_radj,pl_orbeccen,ttv_flag,' +\
        'st_teff,st_rad,st_mass,st_met,st_logg,pl_tranmid,pl_trandur,' +\
        'st_tefferr1,st_raderr1,st_meterr1,st_loggerr1,' +\
        'pl_orbincl,pl_orblper' +\
        '+from+pscomppars&format=csv',
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # EU
        'http://exoplanet.eu/catalog/csv',
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # OEC
        'https://raw.githubusercontent.com/OpenExoplanetCatalogue/' + \
        'oec_tables/master/comma_separated/open_exoplanet_catalogue.txt',
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # ORG
        'http://exoplanets.org/csv-files/exoplanets.csv',
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # TEP
        'https://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv'
    ]
    archive = ['nasa', 'eu', 'oec', 'org', 'tep']
    for i, download_link in enumerate(download_links):
        i = archive[i]
        csv = f'firefly/data/{i}.csv.gz'
        if not os.path.exists(csv):
            print(f'Caching the {i.upper()} Exoplanet Archive.')
            df = read_csv(download_link)
            df.to_csv(csv, index=False)
        two_days_ago = datetime.now() - timedelta(days=2)
        filetime = datetime.fromtimestamp(os.path.getctime(csv))
        if filetime < two_days_ago:
            print(f'{i.upper()} Archive is 2 days old. Updating.')
            df = read_csv(download_link)
            df.to_csv(csv, index=False)
        else:
            pass


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


def _IQR(df, sigma=1):
    # Only keep IQR of data
    if sigma==1:
        Q1 = df.quantile(0.45)
        Q3 = df.quantile(0.55)
    if sigma==2:
        Q1 = df.quantile(0.35)
        Q3 = df.quantile(0.65)
    if sigma==3:
        Q1 = df.quantile(0.25)
        Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    trueList = ~((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR)))
    return df[trueList]
  
    
def priors(exoplanet, archive='eu', save=False, user=True):
    '''
    Generates priors from 4 exoplanet archives, nasa, eu, oec and exoplanets.org

    Parameters
    ----------
    exoplanet : str, 'wasp43b'
        The exoplanet to create a prior object for.
    archive : str, {'eu', 'nasa', 'org', 'all'}
        The archive to generate priors from. All takes the IQR of all
        archives (including OEC) and then the mean.
        The default is 'eu'.
    save : bool, optional
        Internal use only. The default is False.
    user : bool, optional
        Internal use only. The default is True.

    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download NASA archive
    _download_archive()
    if user==True:
        _load_csv()
        highest, ratios = _search_all(exoplanet)
        exoplanet = highest[0]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    archive_list = [exo_nasa, exo_eu, exo_org, exo_oec]
    exo_archive = concat(archive_list).set_index('pl_name')
    # if archive =='all':
    #     exo_archive = concat(archive_list)
    #     exo_archive = exo_archive[exo_archive['pl_name'].notna()]
    #     exo_archive['pl_name'] = \
    #         Categorical(exo_archive['pl_name'],
    #         ordered=True,
    #         categories=natsorted(exo_archive['pl_name'].unique()))
    #     exo_archive = exo_archive.sort_values('pl_name').set_index('pl_name')
    #     exo = exo_archive.mean(axis=0, level=0)
    #     exo.to_csv('firefly/data/all.csv.gz')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    try:
        df = DataFrame(exo_archive).loc[[exoplanet]]
    except KeyError:
        sys.exit('The chosen target is either spelt incorrectly, or does not '
                 'exist in the EU archive.')
    try:
        tic = df['tic_id'] .drop_duplicates() .dropna()[0]
    except IndexError:
        tic = 'N/A'
    # Choose Archive
    if archive=='nasa':
        df = df.iloc[[0]]
    elif archive=='eu':
        df = df.iloc[[1]]
    # elif archive=='oec':
    #     df = df.iloc[[1]]
    elif archive=='org':
        df = df.iloc[[2]]
    # Turn into a series
    s = df.iloc[0]
    t0 = s.loc['pl_tranmid']
    if archive=='all':
        t0 = df['pl_tranmid'].max()
        df = _IQR(df, sigma=3)
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
    rserr = s.loc['st_raderr1']
    z = s.loc['st_met']
    zerr =  s.loc['st_meterr1']
    ms = s.loc['st_mass']
    T = s.loc['st_teff']
    Terr =  s.loc['st_tefferr1']
    t14 = s.loc['pl_trandur'] * 60
    logg = s.loc['st_logg']
    loggerr = s.loc['st_loggerr1']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    host_T = (T, Terr)
    host_z = (z, zerr)
    host_r = (rs, rserr)
    host_logg = (logg, loggerr)
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
    if np.isnan(ecc):
        ecc = 0
    # t14
    if np.isnan(t14):
        t14 = estimate_t14(rp, rs, a, P)
    # logg
    if np.isnan(logg):
        logg, err_logg = calculate_logg((ms, ms * 1e-2), (rs, rs * 1e-2))
        host_logg = (logg, err_logg)
    # z
    if (np.isnan(zerr) and z!=0):
        host_z = (z, 0.1)
    if (np.isnan(z) or z==0):
        host_z = (0, 0.1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Exoplanet Priors to TransitFit
    radius_const = 0.1027626851
    cols = [['P', 'gaussian', P, P * 1e-4, ''],
            ['t0', 'gaussian', t0, 7e-3, ''],
            ['a', 'gaussian', a, a * 1e-2, ''],
            ['inc', 'gaussian', i, i * 1e-2, ''],
            ['w', ['fixed' if (w==90 or w==0) else 'gaussian'][0], w,
             ['' if (w==90 or w==0) else np.abs(w * 0.5)][0], ''],
            ['ecc', ['fixed' if ecc==0 else 'gaussian'][0], ecc,
             ['' if ecc==0 else ecc * 0.5][0], ''],
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
             ['' if (w==90 or w==0) else np.abs(w * 0.5)][0], ''],
            ['ecc', ['fixed' if ecc==0 else 'gaussian'][0], ecc,
             ['' if ecc==0 else ecc * 0.5][0], ''],
            ['rp', 'uniform',
             0.9 * radius_const * rp / rs,
             1.1 * radius_const * rp / rs, 0],
            ['host_T', 'fixed', host_T[0], host_T[1], ''],
            ['host_z', 'fixed', host_z[0], host_z[1], ''],
            ['host_r', 'fixed', host_r[0], host_r[1], ''],
            ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input A', 'Input B', 'Filter'])
    is_nan = repack.isnull().values.any()
    if is_nan:
        raise NaNError(f'Skipping {exoplanet} due to missing prior data.')
    if archive=='all':
        print(f'\nPriors generated from the NASA, EU, OEC and ORG Archives for'
              f' {exoplanet} ({tic}).\n')
    elif archive not in ['nasa', 'eu', 'org']:
        archive = 'NASA'
    else:
        print(f'\nPriors generated from the {archive.upper()} Archive for'
              f' {exoplanet} ({tic}).\n')
    print(tabulate(repack, tablefmt='psql', showindex=False, headers='keys'))
    if user==False:
        return host_T, host_z, host_r, host_logg, t0, P, t14, repack


def tess(archive='eu', survey=None):
    '''
    Currently there are 420 tess targets with full prior and lightcurve sets.

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
    tess = f'{here}/data/Filters/{archive}_tess_viable.csv'
    tess_ttv = f'{here}/data/Filters/tess_ttv_viable.csv'
    targets = read_csv(tess)['Exoplanet']
    ttv_targets = read_csv(tess_ttv)['Exoplanet']
    if survey != None:
        targets = [s for s in targets if survey in s]
        ttv_targets = [s for s in ttv_targets if survey in s]
    all_targets = natsorted(targets)
    ttv_targets = natsorted(ttv_targets)
    return targets, all_targets, ttv_targets


def gen_tess(archive='eu'):
    _download_archive()
    _load_csv()
    exo_list = exo_nasa[['pl_name', 'tic_id']] \
              .dropna() .drop_duplicates('pl_name') \
              .drop(['tic_id'], axis=1) .values .tolist()
    exo_list = [j for i in exo_list for j in i]
    exo_list = natsorted(exo_list)
    viable = []
    products = []
    period = []
    epochs = []
    tic = []
    for i, exoplanet in enumerate(exo_list):
        try:
            host_T, host_z, host_r, host_logg, t0, P, t14, repack = \
                            priors(exoplanet, archive, user=False)
            lc_links, tic_id = _lc(exoplanet)
            if not len(lc_links)==0:
                epoch = ceil((0.8 * 27.4 / P) * len(lc_links))
                viable.append(exoplanet)
                tic.append(tic_id)
                products.append(len(lc_links))
                period.append(P)
                epochs.append(epoch)
            else:
                pass
        except Exception:
            pass
    data = {'Exoplanet':viable, 'TIC ID':tic, 'Products':products, 'Period':period,
            'Epochs':epochs}
    df = DataFrame(data)
    df['Exoplanet'] = \
            Categorical(df['Exoplanet'],
            ordered=True,
            categories=natsorted(df['Exoplanet'].unique()))
    df = df.sort_values('Exoplanet')
    here = os.path.dirname(os.path.abspath(__file__))
    df.to_csv(f'{here}/data/Filters/{archive}_tess_viable.csv', index=False)



def gen_tess_ttv():
    _download_archive()
    _load_csv()
    ttv_list = exo_nasa[['pl_name', 'tic_id', 'ttv_flag']]
    ttv_list = ttv_list[ttv_list!=0] . dropna() \
              .drop_duplicates('pl_name') \
              .drop(['tic_id', 'ttv_flag'], axis=1) .values .tolist()
    ttv_list = [j for i in ttv_list for j in i]
    ttv_list = natsorted(ttv_list)
    viable_ttv = []
    products_ttv = []
    period_ttv = []
    epochs_ttv = []
    tic = []
    for i, exoplanet in enumerate(ttv_list):
        try:
            host_T, host_z, host_r, host_logg, t0, P, t14, repack = \
                            priors(exoplanet, archive, user=False)
            lc_links, tic_id = _lc(exoplanet)
            if not len(lc_links)==0:
                epoch = ceil((0.8 * 27.4 / P) * len(lc_links))
                viable_ttv.append(exoplanet)
                tic.append(tic_id)
                products_ttv.append(len(lc_links))
                period_ttv.append(P)
                epochs_ttv.append(epoch)
            else:
                pass
        except Exception:
            pass
    data = {'Exoplanet':viable_ttv, 'Products':products_ttv,
            'Period':period_ttv, 'Epochs':epochs_ttv}
    df = DataFrame(data)
    df['Exoplanet'] = \
            Categorical(df['Exoplanet'],
            ordered=True,
            categories=natsorted(df['Exoplanet'].unique()))
    df = df.sort_values('Exoplanet')
    here = os.path.dirname(os.path.abspath(__file__))
    df.to_csv(f'{here}/data/Filters/tess_ttv_viable.csv', index=False)
    
