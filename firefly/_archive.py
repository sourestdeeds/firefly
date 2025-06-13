#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exoplanet archive tools.

@author: Stephen Charles
"""

from rich import print
from transitfit import calculate_logg
from datetime import datetime, timedelta
from tabulate import tabulate
from pandas import DataFrame, read_csv, Categorical, concat
from astroquery.mast import Observations as obs
from thefuzz import process
from natsort import natsorted
from tqdm import tqdm
import numpy as np
import time
import sys
import os
import threading

tqdm.pandas(desc="Progress")


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
    download_archive_threaded_execute()
    nasa_csv = 'firefly/data/nasa.csv.gz'
    eu_csv = 'firefly/data/eu.csv.gz'
    oec_csv = 'firefly/data/oec.csv.gz'
    org_csv = 'firefly/data/org.csv.gz'
    here = os.path.dirname(os.path.abspath(__file__))
    spearnet_csv = f'{here}/data/spear.csv'
    spearnet_csv_uncoup = f'{here}/data/spear_uncoupled.csv'
    global exo_nasa, exo_eu, exo_oec, exo_org, exo_spearnet, exo_spearnet_uncoup
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # NASA
    exo_nasa = read_csv(nasa_csv)
    cols = exo_nasa.columns.tolist()
    blank = DataFrame(columns=cols)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # SPEARNET
    try:
        exo_spearnet = read_csv(spearnet_csv)
        exo_spearnet = concat([exo_spearnet, blank], ignore_index=True)
        exo_spearnet['archive'] = 'SPEARNET COUPLED'
    except FileNotFoundError:
        print(f"SPEARNET coupled file not found: {spearnet_csv}")
        exo_spearnet = DataFrame(columns=cols + ['archive'])
        exo_spearnet['archive'] = 'SPEARNET COUPLED'

    try:
        exo_spearnet_uncoup = read_csv(spearnet_csv_uncoup)
        exo_spearnet_uncoup = concat([exo_spearnet_uncoup, blank], ignore_index=True)
        exo_spearnet_uncoup['archive'] = 'SPEARNET UNCOUPLED'
    except FileNotFoundError:
        print(f"SPEARNET uncoupled file not found: {spearnet_csv_uncoup}")
        exo_spearnet_uncoup = DataFrame(columns=cols + ['archive'])
        exo_spearnet_uncoup['archive'] = 'SPEARNET UNCOUPLED'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # EU
    try:
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
    except ValueError:
        # Handle case where column names have changed
        try:
            col_subset_eu = ['name', 'orbital_period', 'semi_major_axis', 'radius',
                             'eccentricity', 'inclination', 'tzero_tr', 'star_teff',
                             'star_teff_error_max', 'star_radius',
                             'star_radius_error_max', 'star_mass',
                             'star_metallicity', 'star_metallicity_error_max', 'omega']
            exo_eu = read_csv(eu_csv, usecols=col_subset_eu)
            exo_eu.columns = ['pl_name', 'pl_radj', 'pl_orbper', 'pl_orbsmax',
                              'pl_orbeccen', 'pl_orbincl', 'pl_orblper', 'pl_tranmid', 'st_met',
                              'st_meterr1', 'st_mass','st_rad',
                              'st_raderr1', 'st_teff', 'st_tefferr1']
        except ValueError:
            # If specific columns still don't work, read all columns and try to map them
            exo_eu_full = read_csv(eu_csv)
            print(f"EU CSV columns: {list(exo_eu_full.columns)}")
            # Create a minimal dataframe with the required structure
            exo_eu = DataFrame(columns=['pl_name', 'pl_radj', 'pl_orbper', 'pl_orbsmax',
                                        'pl_orbeccen', 'pl_orbincl', 'pl_orblper', 'pl_tranmid', 'st_met',
                                        'st_meterr1', 'st_mass','st_rad',
                                        'st_raderr1', 'st_teff', 'st_tefferr1'])
    exo_eu = concat([exo_eu, blank], ignore_index=True)
    exo_eu['archive'] = 'EU'
    exo_nasa['archive'] = 'NASA'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # OEC
    try:
        col_subset_oec = ['name', 'radius', 'period',
                         'semimajoraxis', 'eccentricity', 'periastron', 'inclination',
                         'hoststar_mass', 'hoststar_radius',
                         'hoststar_metallicity', 'hoststar_temperature']
        exo_oec = read_csv(oec_csv, usecols=col_subset_oec)
        exo_oec.columns = ['pl_name', 'pl_radj', 'pl_orbper', 'pl_orbsmax',
                          'pl_orbeccen', 'pl_orblper','pl_orbincl', 'st_mass',
                          'st_rad', 'st_met', 'st_teff']
    except ValueError:
        # Handle case where column names have changed or don't exist
        exo_oec_full = read_csv(oec_csv)
        print(f"OEC CSV columns: {list(exo_oec_full.columns)}")
        # Create a minimal dataframe with the required structure
        exo_oec = DataFrame(columns=['pl_name', 'pl_radj', 'pl_orbper', 'pl_orbsmax',
                                     'pl_orbeccen', 'pl_orblper','pl_orbincl', 'st_mass',
                                     'st_rad', 'st_met', 'st_teff'])
    exo_oec = concat([exo_oec, blank], ignore_index=True)
    exo_oec['archive'] = 'OEC'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # ORG
    try:
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
    except ValueError:
        # Handle case where column names have changed or don't exist
        exo_org_full = read_csv(org_csv)
        print(f"ORG CSV columns: {list(exo_org_full.columns)}")
        # Create a minimal dataframe with the required structure
        exo_org = DataFrame(columns=['pl_orbeccen', 'st_met', 'st_meterr1', 'pl_orbincl',
                                     'st_logg', 'st_loggerr1','st_mass', 'pl_name',
                                     'pl_orbper', 'pl_radj', 'st_rad', 'st_raderr1', 'pl_orbsmax',
                                     'st_teff', 'st_tefferr1', 'pl_tranmid'])
    exo_org = concat([exo_org, blank], ignore_index=True)
    exo_org['archive'] = 'ORG'
    

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
    exo_archive = concat(archive_list)
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

def _pl(tic_id):
    '''
    Returns a dataframe of all planet names and tic ID's.

    Returns
    -------
    Planet name corresponding to tic_id
    '''
    df = exo_nasa[['pl_name','tic_id']].dropna() \
        .sort_values('pl_name')
    df['tic_id'] = df['tic_id'].str[3:]
    df['tic_id'] = df['tic_id'].astype(int)
    df = df.set_index('tic_id')
    try:
        target = df.loc[[tic_id]]
    except KeyError:
        target = ''
    try:
        exoplanet = target['pl_name'] .values .tolist()[0]
    except Exception:
        exoplanet = ''
    return exoplanet


def _lc(exoplanet, mast, fast=False):
    '''
    # https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html
    out = DataFrame()
    for i in range(1,32):
        a = read_csv(f'TESS_curl/tesscurl_sector_{str(i)}_lc.sh')
        a['links'] = a['#!/bin/sh'].str.split().str[-1]
        a = a.drop(['#!/bin/sh'], axis=1)
        out = concat([out, a], ignore_index=True)
    out.to_csv('MAST_lc.csv.gz', index=False)
    '''
    tic_id = _tic(exoplanet).replace('TIC ', '')
    lc_links = mast[mast['links'].str.contains(tic_id)] .values .tolist()
    lc_list = [i for j in lc_links for i in j]
    if fast==True:
        lc_test = [int(i[-35:-20]) for j in lc_links for i in j]
    else:
        lc_test = [int(i[-30:-15]) for j in lc_links for i in j]
    lc_links = []
    for i in range(len(lc_test)):
        if int(tic_id) == lc_test[i]:
            lc_links.append(lc_list[i])
    return lc_links, tic_id


def tess_candidate():
    '''
    Finds all candidate tic_ids and groups by amount of sector observations.

    '''
    _load_csv()
    here = os.path.dirname(os.path.abspath(__file__))
    mast_csv = f'{here}/data/Search/TESS_lc.csv.xz'
    mast = read_csv(mast_csv)
    for i in range(13,21):
        count = i
        df = mast.sort_values('links').reset_index()
        df = df.groupby(df.links.str[-30:-15].astype(int)).size().reset_index(name='Products')
        df = df.loc[df['Products']==count]
        df = df.rename(columns={'links':'TIC ID'})
        df['Exoplanet'] = df.apply(lambda row: _pl(row['TIC ID']), axis=1)
        here = os.path.dirname(os.path.abspath(__file__))
        os.makedirs(f'{here}/data/Candidates', exist_ok=True)
        df.to_csv(f'{here}/data/Candidates/{count}_sector_candidates.csv', index=False)
    

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
        'pl_name,tic_id,pl_orbper,pl_orbsmax,pl_radj,pl_bmasse,pl_bmassj,pl_orbeccen,ttv_flag,' +\
        'st_teff,st_rad,st_mass,st_met,st_logg,pl_tranmid,pl_trandur,' +\
        'st_tefferr1,st_raderr1,st_meterr1,st_loggerr1,' +\
        'pl_orbincl,pl_orblper,ra,dec,glat,glon,sy_dist,sy_plx,sy_tmag,' +\
        'pl_orbpererr1,pl_tranmiderr1,pl_orbsmaxerr1,pl_radjerr1,pl_orbinclerr1,' +\
        'pl_orblpererr1,pl_orbeccenerr1'
        '+from+pscomppars&format=csv',
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # NASA Kepler Names
        'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/' +\
        'nstedAPI/nph-nstedAPI?' +\
        'table=k2pandc',
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
    archive = ['nasa', 'eu', 'kep', 'oec', 'org', 'tep']
    for i, download_link in enumerate(download_links):
        try:
            i = archive[i]
            csv = f'firefly/data/{i}.csv.gz'
            if not os.path.exists(csv):
                print(f'Caching the {i.upper()} Exoplanet Archive.')
                df = read_csv(download_link)
                df.to_csv(csv, index=False)
            seven_days_ago = datetime.now() - timedelta(days=7)
            filetime = datetime.fromtimestamp(os.path.getctime(csv))
            if filetime < seven_days_ago:
                print(f'{i.upper()} Archive is 7 days old. Updating.')
                df = read_csv(download_link)
                df.to_csv(csv, index=False)
            else:
                pass
        except Exception:
            print(f'Could not contact the server for the {i.upper()} archive.')
            pass


def _download_archive_threaded(download_link, i):
    os.makedirs('firefly/data', exist_ok=True)
    try:
        csv = f'firefly/data/{i}.csv.gz'
        if not os.path.exists(csv):
            print(f'Caching the {i.upper()} Exoplanet Archive.')
            df = read_csv(download_link)
            df.to_csv(csv, index=False)
        seven_days_ago = datetime.now() - timedelta(days=7)
        filetime = datetime.fromtimestamp(os.path.getctime(csv))
        if filetime < seven_days_ago:
            print(f'{i.upper()} Archive is 7 days old. Updating.')
            df = read_csv(download_link)
            df.to_csv(csv, index=False)
        else:
            pass
    except Exception:
        print(f'Could not contact the server for the {i.upper()} archive.')
        pass

def download_archive_threaded_execute():
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
        'pl_name,tic_id,pl_orbper,pl_orbsmax,pl_radj,pl_bmasse,pl_bmassj,pl_orbeccen,ttv_flag,' +\
        'st_teff,st_rad,st_mass,st_met,st_logg,pl_tranmid,pl_trandur,' +\
        'st_tefferr1,st_raderr1,st_meterr1,st_loggerr1,' +\
        'pl_orbincl,pl_orblper,ra,dec,glat,glon,sy_dist,sy_plx,sy_tmag,' +\
        'pl_orbpererr1,pl_tranmiderr1,pl_orbsmaxerr1,pl_radjerr1,pl_orbinclerr1,' +\
        'pl_orblpererr1,pl_orbeccenerr1'
        '+from+pscomppars&format=csv',
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # NASA Kepler Names
        'https://exoplanetarchive.ipac.caltech.edu/' +\
        'TAP/sync?query=select+*+from+k2pandc&format=csv',
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
    archive = ['nasa', 'kep', 'eu', 'oec', 'org', 'tep']
    threads = []
    for i in range(5):
        j = archive[i]
        download_link = download_links[i]
        t = threading.Thread(target=_download_archive_threaded, args=(download_link, j))
        t.daemon = True
        threads.append(t)
    for i in range(5):
        threads[i].start()
    for i in range(5):    
        threads[i].join()


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
  
    
def priors(exoplanet, archive='eu', save=False, user=True, auto=True, fit_ttv=False):
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
    #_download_archive()
    if user==True:
        download_archive_threaded_execute()
        _load_csv()
        highest, ratios = _search_all(exoplanet)
        exoplanet = highest[0]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if fit_ttv==True:
        archive = 'spearnet'
    archive_list = [exo_nasa, exo_eu, exo_org, exo_oec, exo_spearnet, exo_spearnet_uncoup]
    if archive=='nasa':
        exo_archive = archive_list[0].set_index('pl_name')
    elif archive=='eu':
        exo_archive = archive_list[1].set_index('pl_name')
    elif archive=='org':
        exo_archive = archive_list[2].set_index('pl_name')
    elif archive=='spearnet':
        exo_archive = archive_list[4].set_index('pl_name')
    elif archive=='all':
        exo_archive = concat(archive_list).set_index('pl_name')
    if archive =='all':
        exo_archive = concat(archive_list)
        exo_archive = exo_archive[exo_archive['pl_name'].notna()]
        exo_archive['pl_name'] = \
            Categorical(exo_archive['pl_name'],
            ordered=True,
            categories=natsorted(exo_archive['pl_name'].unique()))
        exo_archive = exo_archive.sort_values('pl_name').set_index('pl_name') \
                      .drop(['Archive'], axis=1)
        # exo = exo_archive.mean(axis=0, level=0)
        exo_archive.to_csv('firefly/data/all.csv.gz')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    try:
        choice = exo_archive.loc[[exoplanet]]
    except KeyError:
        sys.exit('The chosen target is either spelt incorrectly, or does not '
                 f'exist in the {archive.upper()} archive.')
    try:
        # tic = df['tic_id'] .drop_duplicates() .dropna()[0]
        tic = _tic(exoplanet)
    except IndexError:
        tic = 'N/A'
    # Choose Archive
    if archive=='spearnet':
        exo_limb = archive_list[0].set_index('pl_name')
        choice_limb = exo_limb.loc[[exoplanet]]
        s_limb = choice_limb.mean()
    # Turn into a series
    s = choice.iloc[0]
    t0 = s.loc['pl_tranmid']
    if archive=='all':
        t0 = choice['pl_tranmid'].max()
        df = _IQR(choice, sigma=3)
        s = df.mean()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Values for calculation
    P = s.loc['pl_orbper']
    a = s.loc['pl_orbsmax']
    i = s.loc['pl_orbincl']
    w = s.loc['pl_orblper']
    ecc = s.loc['pl_orbeccen']
    rp = s.loc['pl_radj']
    # rserr = s.loc['st_raderr1']
    rs = s.loc['st_rad']
    z = s.loc['st_met']
    zerr =  s.loc['st_meterr1']
    ms = s.loc['st_mass']
    T = s.loc['st_teff']
    # Terr =  s.loc['st_tefferr1']
    t14 = s.loc['pl_trandur'] * 60
    logg = s.loc['st_logg']
    # loggerr = s.loc['st_loggerr1']
    if archive=='nasa':
        ra = s.loc['ra']
        dec = s.loc['dec']
        dist = s.loc['sy_dist']
        t_mag = s.loc['sy_tmag']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    host_T = (T, T * 5e-2)
    host_z = (z, zerr)
    host_r = (rs, rs * 5e-2)
    host_logg = (logg, logg * 5e-2)
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
    if archive=='spearnet':
        Perr = s.loc['pl_orbpererr1']
        aerr = s.loc['pl_orbsmaxerr1']
        ierr = s.loc['pl_orbinclerr1']
        werr = s.loc['pl_orblpererr1']
        eccerr = s.loc['pl_orbeccenerr1']
        t0err = s.loc['pl_tranmiderr1']
        rperr = s.loc['pl_radjerr1']
        rs = s_limb.loc['st_rad']
        z = s_limb.loc['st_met']
        zerr =  s_limb.loc['st_meterr1']
        ms = s_limb.loc['st_mass']
        T = s_limb.loc['st_teff']
        logg = s_limb.loc['st_logg']
        host_T = (T, T * 5e-2)
        host_z = (z, zerr)
        host_r = (rs, rs * 5e-2)
        host_logg = (logg, logg * 5e-2)
        t14 = estimate_t14(rp, rs, a, P)
    if (np.isnan(w) or w==0):
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
        host_z = (z, 0.05)
    if (np.isnan(z) or z==0):
        host_z = (0, 0.05)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Exoplanet Priors to TransitFit
    radius_const = 0.1027626851
    cols = [['P', ['fixed' if archive=='spearnet' else 'gaussian'][0], P,
                ['' if archive=='spearnet' else P * 1e-3][0], ''],
            ['t0', 'gaussian', t0, [7e-3 if archive=='spearnet' else 7e-3][0], ''],
            ['a', 'gaussian', a, [a * 2e-1 if archive=='spearnet' else a * 2e-1][0], ''],
            ['inc', 'gaussian', i, [i * 2e-1 if archive=='spearnet' else i * 2e-1][0], ''],
            ['w', ['fixed' if (w==90 or ecc==0) else 'gaussian'][0],
             [90 if ecc==0 else w][0],
             ['' if (w==90 or ecc==0)
              else (360-np.abs(w))/2 if 180 < w < 270
              else (360-np.abs(w)) if 270 <= w < 360
              else werr if archive=='spearnet'
              else w * 0.5][0], ''],
            ['ecc', ['fixed' if (ecc==0 or w==90)
                     else 'gaussian'][0], [0 if w==90 else ecc][0],
             ['' if (ecc==0 or w==90)
              else ecc * 0.5 if archive=='spearnet'
              else ecc * 0.5][0], ''],
            ['rp', ['uniform' if archive=='spearnet' else 'uniform'][0],
             [0.5*rp if archive=='spearnet' else (0.5 * radius_const * rp / rs)][0],
             [1.5*rp if archive=='spearnet' else (1.5 * radius_const * rp / rs)][0], 0]]
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
    cols = [['P', ['fixed' if archive=='spearnet' else 'gaussian'][0], P,
                ['' if archive=='spearnet' else P * 1e-3][0], ''],
            ['t0', 'gaussian', t0, [7e-3 if archive=='spearnet' else 7e-3][0], ''],
            ['a', 'gaussian', a, [a * 2e-1 if archive=='spearnet' else a * 2e-1][0], ''],
            ['inc', 'gaussian', i, [i * 2e-1 if archive=='spearnet' else i * 2e-1][0], ''],
            ['w', ['fixed' if (w==90 or ecc==0) else 'gaussian'][0],
             [90 if ecc==0 else w][0],
             ['' if (w==90 or ecc==0)
              else (360-np.abs(w))/2 if 180 < w < 270
              else (360-np.abs(w)) if 270 <= w < 360
              else werr if archive=='spearnet'
              else w * 0.5][0], ''],
            ['ecc', ['fixed' if (ecc==0 or w==90)
                     else 'gaussian'][0], [0 if w==90 else ecc][0],
             ['' if (ecc==0 or w==90)
              else ecc * 0.5 if archive=='spearnet'
              else ecc * 0.5][0], ''],
            ['rp', ['uniform' if archive=='spearnet' else 'uniform'][0],
             [0.5*rp if archive=='spearnet' else (0.5 * radius_const * rp / rs)][0],
             [1.5*rp if archive=='spearnet' else (1.5 * radius_const * rp / rs)][0], 0],
            ['t14', '', t14, '', ''],
            ['host_T', '', host_T[0], host_T[1], ''],
            ['host_z', '', host_z[0], host_z[1], ''],
            ['host_r', '', host_r[0], host_r[1], ''],
            ['host_logg', '', host_logg[0], host_logg[1], '']]
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input A', 'Input B', 'Filter'])
    is_nan = repack.isnull().values.any()
    #if (is_nan and user==False and auto==True):
     #   raise NaNError(f'Skipping {exoplanet} due to missing prior data.')
    if archive=='all':
        print(f'\nPriors generated from the NASA, EU, OEC and ORG Archives for'
              f' {exoplanet} ({tic}).\n')
    elif archive not in ['nasa', 'eu', 'org', 'spearnet']:
        archive = 'NASA'
    else:
        print(f'\nPriors generated from the {archive.upper()} Archive for'
              f' {exoplanet} ({tic}).\n')
    print(tabulate(repack, tablefmt='psql', showindex=False, headers='keys'))
    if (archive=='nasa' and user==False):
        return host_T, host_z, host_r, host_logg, t0, P, t14, repack, ra, dec, dist, t_mag
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
    tess = f'{here}/data/Targets/{archive}_tess_viable.csv'
    tess_ttv = f'{here}/data/Targets/tess_ttv_viable.csv'
    targets = read_csv(tess)['Exoplanet']
    ttv_targets = read_csv(tess_ttv)['Exoplanet']
    if survey != None:
        targets = [s for s in targets if survey in s]
        ttv_targets = [s for s in ttv_targets if survey in s]
    all_targets = natsorted(targets)
    ttv_targets = natsorted(ttv_targets)
    return targets, all_targets, ttv_targets


def gen_tess(archive='nasa', cadence=120):
    '''Generates all TESS targets which have a fitsfile.'''

    def count_products(tic_id):
        '''Counts the amount of TESS fitsfiles on MAST.'''
        time.sleep(0.01)
        search = obs.query_criteria(dataproduct_type=['timeseries'],
                                    project='TESS',
                                    provenance_name=['SPOC'],
                                    t_exptime=cadence,
                                    target_name=tic_id).to_pandas()
        search = search[search['t_exptime']==cadence]
        search = search[~search.dataURL.str.endswith('_dvt.fits')]
        if len(search) is None:
            return 0
        return len(search)

    _load_csv()
    if archive=='nasa':
        df = exo_nasa
    elif archive=='eu':
        df = exo_eu
    elif archive=='nasa_full':
        download_link =  \
        'https://exoplanetarchive.ipac.caltech.edu/' +\
        'TAP/sync?query=select+*+from+ps&format=csv'
        print('Downloading the full NASA archive.')
        df = read_csv(download_link)
        print('Download complete.')
    print('Generating product list and transit counts.')
    print('The process should take no longer than 15 minutes.')
    df['tic_id'] = df['tic_id'].str.replace('TIC ', '')
    df[f'Products ({cadence} Cadence)'] = df['tic_id'].progress_apply(count_products)
    df['Transits'] = np.ceil((0.8 * 27.4 / df['pl_orbper']) * df[f'Products ({cadence} Cadence)'])
    df['Transits'] = df['Transits'].fillna(-1).astype(int)
    df['pl_name'] = \
            Categorical(df['pl_name'],
            ordered=True,
            categories=natsorted(df['pl_name'].unique()))
    df = df.sort_values('pl_name')
    here = os.path.dirname(os.path.abspath(__file__))
    df.to_csv(f'{here}/data/Targets/ML_nasa_tess_viable.csv.xz', index=False)
    df.to_csv('firefly/data/ML_nasa_tess_viable.csv.xz', index=False)
    print('Operation completed.')

