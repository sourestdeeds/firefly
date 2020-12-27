#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data retrievers.

@author: Steven Charles-Mindoza
"""


from ._search import _fuzzy_search
from ._archive import _eu, _download_nasa, _nasa

try:
    from lightkurve import search_lightcurve
except:
    try:
        from lightkurve import search_lightcurvefile as search_lightcurve
    except:
        raise
from tabulate import tabulate
from shutil import rmtree
from pandas import read_csv
import sys
import os



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
    nasa_csv = 'firefly/data/nasa.csv'
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
    nasa_csv = 'firefly/data/nasa.csv'
    df = read_csv(nasa_csv, usecols=['pl_name',
                                     'tic_id',
                                     'soltype',
                                     'rowupdate',
                                     'disc_facility']) .sort_values('pl_name') \
                                    .drop_duplicates('pl_name', keep='last') 
                                    #.sort_values('pl_name')
    tic = 'firefly/data/tic.csv'
    df.to_csv(tic, index=False, header=True)
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
    highest, ratios = _fuzzy_search(target, archive=archive)
    exoplanet = highest[0]
    lc = search_lightcurve(exoplanet, mission='TESS')
    if len(lc) == 0:
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    lc = lc .table .to_pandas()[['observation', 
                                 'productFilename', 'size', 't_exptime']] \
            .rename(columns={'observation':'Observation'}) \
            .rename(columns={'size':'Size'}) \
            .rename(columns={'productFilename':'Product'}) 
    lc = lc[lc.t_exptime != 20].drop(['t_exptime'], axis=1)
    print(tabulate(lc, tablefmt='psql', showindex=False, headers='keys'))
    
