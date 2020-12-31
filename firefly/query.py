#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data retrievers.

@author: Steven Charles-Mindoza
"""


from ._search import _fuzzy_search
from ._archive import _eu, _download_nasa, _nasa

from tabulate import tabulate
from shutil import rmtree
from pandas import read_csv, DataFrame
import random
import sys
import os


def tess_viable(k=10):
    '''
    Currently there are 377 tess targets with full prior and lightcurve sets.

    Parameters
    ----------
    k : int, optional
        Returns a random set of targets. The default is 10.

    Returns
    -------
    targets : str
        Viable tess target list.

    '''
    here = os.path.dirname(os.path.abspath(__file__))
    tess = f'{here}/data/Filters/tess_viable.csv'
    targets = read_csv(tess) 
    targets = targets['Exoplanet'] .values .tolist()
    targets = random.sample(targets, k=k)
    return targets


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
    here = os.path.dirname(os.path.abspath(__file__))
    mast = f'{here}/data/Filters/MAST_lc.csv.gz'
    _ = read_csv(mast)
    tic_id = tic(exoplanet).replace('TIC ', '')
    lc_links = _[_['links'].str.contains(tic_id)] .values .tolist()
    lc_links = [i for j in lc_links for i in j]
    return lc_links, tic_id
    
    

def query(target, archive='eu'):
    '''
    Performs a name search and returns the 5 top matches. 
    Parameters
    ----------
    exoplanet : str
        The exoplanet target to retreive information for.
    archive : str, optional
        The exoplanet archive to use for priors. The default is 'eu'.

    Returns
    -------
    Data printed to console.

    '''
    output = _fuzzy_search(target, archive=archive)
    table = output[1]
    print(tabulate(table, tablefmt='psql', headers=['Exoplanet', '% Match']))


def tess():
    '''
    A list of all exoplanets with a TIC ID, and full prior entries.

    Parameters
    ----------
    None

    Returns
    -------
    Data printed to console.

    '''
    _download_nasa()
    nasa_csv = 'firefly/data/nasa.csv.gz'
    df = read_csv(nasa_csv, usecols=['pl_name',
                                     'tic_id',
                                     'soltype',
                                     'rowupdate',
                                     'disc_facility']) \
                                    .drop_duplicates('pl_name', keep='last') \
                                    .sort_values('pl_name')
    TESS = 'Transiting Exoplanet Survey Satellite (TESS)'
    disc_by_tess = df[df.disc_facility == TESS].drop(['disc_facility'], axis=1)
    print(tabulate(disc_by_tess, tablefmt='psql', showindex=False,
                   headers=['Exoplanet', 'TIC ID', 'Confirmed/Candidate', 'Updated']))
    tess_targets = disc_by_tess .drop(['tic_id', 'soltype', 'rowupdate'], axis=1) \
                 .values .tolist()
    tess_targets = [i for j in tess_targets for i in j]
            
    return tess_targets


def tic(exoplanet):
    '''
    Returns a dataframe of all planet names and tic ID's.

    Parameters
    ----------
    None

    Returns
    -------
    TIC ID of an exoplanet target.

    '''
    _download_nasa()
    nasa_csv = 'firefly/data/nasa.csv.gz'
    df = read_csv(nasa_csv, usecols=['pl_name',
                                     'tic_id',
                                     'soltype',
                                     'rowupdate',
                                     'disc_facility']) .sort_values('pl_name') \
                                    .drop_duplicates('pl_name', keep='last')
    tic_df = df .drop(['soltype', 'rowupdate', 'disc_facility'], axis=1).dropna() \
                 .set_index('pl_name')
    target = tic_df.loc[[exoplanet]]
    tic_id = target['tic_id'] .values .tolist()[0]
    return tic_id


def priors(target, archive='eu'):
    '''
    Performs a query for prior information.

    Parameters
    ----------
    exoplanet : str
        The exoplanet target to retreive information for.
    archive : str, optional
        The exoplanet archive to use for priors. The default is 'eu'.

    Returns
    -------
    Data printed to console.

    '''
    highest, ratios = _fuzzy_search(target, archive=archive)
    exoplanet = highest[0]
    temp = f'firefly/{exoplanet}'
    os.makedirs(temp, exist_ok=True)
    if archive == 'eu':
        _eu(exoplanet)
    elif archive == 'nasa':
        _nasa(exoplanet)
    rmtree(temp)

def mast(target, archive='eu'):
    '''
    Performs a query for data products from MAST.

    Parameters
    ----------
    exoplanet : str
        The exoplanet target to retreive information for.
    archive : str, optional
        The exoplanet archive to use for priors. The default is 'eu'.

    Returns
    -------
    Data printed to console.

    '''
    exo_folder = f'firefly/{target}'
    highest, ratios = _fuzzy_search(target, archive=archive)
    exoplanet = highest[0]
    print(f'\nSearching MAST for {exoplanet}.')
    lc_links, tic_id = _lc(exoplanet)
    if len(lc_links) == 0:
        rmtree(exo_folder)
        print(f'Search result contains no data products for {exoplanet}.')
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    print(f'\nQuery from MAST returned {len(lc_links)} '
          f'data products for {exoplanet} (TIC {tic_id}).\n')
    return lc_links
    
