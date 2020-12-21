from transitfit import calculate_logg
from datetime import datetime, timedelta
from tabulate import tabulate
from pandas import DataFrame, read_csv
from shutil import rmtree
import astropy.units as u
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
        

def _download_eu():
    os.makedirs('firefly/data', exist_ok=True)
    eu_csv = 'firefly/data/eu.csv'
    download_link = 'http://exoplanet.eu/catalog/csv'
    if not os.path.exists(eu_csv):
        print('eu.csv does not exist, downloading.')
        df = read_csv(download_link)
        df.to_csv(eu_csv, index=False)
    ten_days_ago = datetime.now() - timedelta(days=10)
    filetime = datetime.fromtimestamp(os.path.getctime(eu_csv))
    if filetime < ten_days_ago:
        print('eu.csv is 10 days old, updating.')
        df = read_csv(download_link)
        df.to_csv(eu_csv, index=False)
    else:
        pass
    return eu_csv
    

def _download_nasa():
    os.makedirs('firefly/data', exist_ok=True)
    download_link =  \
        'https://exoplanetarchive.ipac.caltech.edu/' +\
        'TAP/sync?query=select+' +\
        'pl_name,pl_orbper,pl_orbpererr1,pl_orbsmax,pl_orbsmaxerr1,' +\
        'pl_radj,pl_radjerr1,pl_orbeccen,pl_orbeccenerr1,disc_facility,' +\
        'st_teff,st_tefferr1,st_rad,st_raderr1,st_mass,st_masserr1,' +\
        'st_met,st_meterr1,pl_tranmid,pl_tranmiderr1,pl_trandur,' +\
        'pl_orbincl,pl_orbinclerr1,pl_orblper,pl_orblpererr1,rowupdate' +\
        '+from+ps&format=csv'
    nasa_csv = 'firefly/data/nasa.csv'
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


def _eu(exoplanet):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download EU archive
    eu_csv = _download_eu()
    col_subset_eu = ['# name', 'orbital_period', 'orbital_period_error_max',
                     'semi_major_axis', 'semi_major_axis_error_max',
                     'radius', 'radius_error_max',
                     'eccentricity', 'eccentricity_error_max',
                     'inclination', 'inclination_error_max',
                     'tzero_tr', 'tzero_tr_error_max',
                     'star_teff', 'star_teff_error_max',
                     'star_radius', 'star_radius_error_max',
                     'star_mass', 'star_mass_error_max',
                     'star_metallicity', 'star_metallicity_error_max']
    exo_archive = read_csv(eu_csv, index_col='# name',
                           usecols=col_subset_eu)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    try:
        df = DataFrame(exo_archive).loc[[exoplanet]]
    except KeyError:
        sys.exit('The chosen target is either spelt incorrectly, or does not '
                 'exist in the EU archive.')
    s = df.mean()
    t0, P, i = s.loc['tzero_tr'], s.loc['orbital_period'], s.loc['inclination']
    rp, rs, a = s.loc['radius'], s.loc['star_radius'], s.loc['semi_major_axis']
    rs_m = (rs * u.R_sun).to(u.m)
    rp_m = (rp * u.R_jup).to(u.m)
    a_m = (a * u.AU).to(u.m)
    t14 = P * (rs_m * np.sin(np.radians(i)) + rp_m) / \
              (np.pi * a_m) * 24 * 60 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    logg, err_logg = calculate_logg((s.loc['star_mass'],
                                     s.loc['star_mass_error_max']),
                                    (s.loc['star_radius'],
                                     s.loc['star_radius_error_max']))
    host_T = (s.loc['star_teff'], s.loc['star_teff_error_max'])
    host_z = (s.loc['star_metallicity'],
              s.loc['star_metallicity_error_max'])
    host_r = (s.loc['star_radius'], s.loc['star_radius_error_max'])
    host_logg = (logg, err_logg)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Exoplanet Priors to TransitFit
    radius_const = 0.1027626851
    cols = [['P', 'gaussian', s.loc['orbital_period'],
             s.loc['orbital_period_error_max'] * 1e5, ''],
            ['t0', 'gaussian', s.loc['tzero_tr'],
             s.loc['tzero_tr_error_max'] * 100, ''],
            ['a', 'gaussian', s.loc['semi_major_axis'],
             s.loc['semi_major_axis_error_max'], ''],
            ['inc', 'gaussian', s.loc['inclination'],
             s.loc['inclination_error_max'], ''],
            ['rp', 'uniform',
             0.5 * radius_const * s.loc['radius'] / s.loc['star_radius'],
             2 * radius_const * s.loc['radius'] / s.loc['star_radius'], 0]]
    priors_csv = DataFrame(cols, columns=['Parameter', 'Distribution',
                                          'Input_A', 'Input_B', 'Filter'])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save the priors
    priors = f'firefly/{exoplanet}/{exoplanet} Priors.csv'
    if not os.path.exists(priors):
        priors_csv.to_csv(priors, index=False, header=True)
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # For printing variables only
    cols = [['P', 'gaussian', s.loc['orbital_period'],
             s.loc['orbital_period_error_max'] * 1e5, ''],
            ['t0', 'gaussian', s.loc['tzero_tr'],
             s.loc['tzero_tr_error_max'] * 100, ''],
            ['a', 'gaussian', s.loc['semi_major_axis'],
             s.loc['semi_major_axis_error_max'], ''],
            ['inc', 'gaussian', s.loc['inclination'],
             s.loc['inclination_error_max'], ''],
            ['rp', 'uniform',
             0.5 * radius_const * s.loc['radius'] / s.loc['star_radius'],
             2 * radius_const * s.loc['radius'] / s.loc['star_radius'], 0],
            ['host_T', 'fixed', host_T[0], host_T[1], ''],
            ['host_z', 'fixed', host_z[0], host_z[1], ''],
            ['host_r', 'fixed', host_r[0], host_r[1], ''],
            ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input_A', 'Input_B', 'Filter'])
    nan = repack.isnull().values.any()
    print(f'\nPriors generated from the EU Archive for {exoplanet}.\n')
    print(tabulate(repack, tablefmt='psql', showindex=False, headers='keys'))
    return host_T, host_z, host_r, host_logg, t0, P, t14, nan


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
    df.to_csv('firefly/data/nasa_full.csv', index=False)
        
    
def _nasa(exoplanet):
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
    s = df.mean()
    t0, P, t14 = s.loc['pl_tranmid'], s.loc['pl_orbper'], \
        s.loc['pl_trandur'] * 60
    print(t14)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    logg, err_logg = calculate_logg((s.loc['st_mass'],
                                     s.loc['st_masserr1']),
                                    (s.loc['st_rad'],
                                     s.loc['st_raderr1']))
    host_T = (s.loc['st_teff'], s.loc['st_tefferr1'])
    host_z = (s.loc['st_met'], s.loc['st_meterr1'])
    host_r = (s.loc['st_rad'], s.loc['st_raderr1'])
    host_logg = (logg, err_logg)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Exoplanet Priors to TransitFit
    radius_const = 0.1027626851
    cols = [['P', 'gaussian', s.loc['pl_orbper'],
             s.loc['pl_orbpererr1'] * 1e5, ''],
            ['t0', 'gaussian', s.loc['pl_tranmid'],
             s.loc['pl_tranmiderr1'] * 100, ''],
            ['a', 'gaussian', s.loc['pl_orbsmax'],
             s.loc['pl_orbsmaxerr1'], ''],
            ['inc', 'gaussian', s.loc['pl_orbincl'],
             s.loc['pl_orbinclerr1'], ''],
            ['rp', 'uniform',
             0.5 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'],
             2 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'], 0]]
    priors_csv = DataFrame(cols, columns=['Parameter', 'Distribution',
                                          'Input_A', 'Input_B', 'Filter'])
    priors_csv = DataFrame(cols, columns=['Parameter', 'Distribution',
                                          'Input_A', 'Input_B', 'Filter'])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save the priors
    priors = f'firefly/{exoplanet}/{exoplanet} Priors.csv'
    if not os.path.exists(priors):
        priors_csv.to_csv(priors, index=False, header=True)
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # For printing variables only
    cols = [['P', 'gaussian', s.loc['pl_orbper'],
             s.loc['pl_orbpererr1'] * 1e5, ''],
            ['t0', 'gaussian', s.loc['pl_tranmid'],
             s.loc['pl_tranmiderr1'] * 100, ''],
            ['a', 'gaussian', s.loc['pl_orbsmax'],
             s.loc['pl_orbsmaxerr1'], ''],
            ['inc', 'gaussian', s.loc['pl_orbincl'],
             s.loc['pl_orbinclerr1'], ''],
            ['rp', 'uniform',
             0.5 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'],
             2 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'], 0],
            ['host_T', 'fixed', host_T[0], host_T[1], ''],
            ['host_z', 'fixed', host_z[0], host_z[1], ''],
            ['host_r', 'fixed', host_r[0], host_r[1], ''],
            ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input_A', 'Input_B', 'Filter'])
    nan = repack.isnull().values.any()
    print(f'\nPriors generated from the NASA Archive for {exoplanet}.\n')
    print(tabulate(repack, tablefmt='psql', showindex=False, headers='keys'))
    return host_T, host_z, host_r, host_logg, t0, P, t14, nan


def _check_nan(exoplanet, archive='eu', printing=False):
    '''
    Checks if the priors generated have values for each entry.

    Parameters
    ----------
    exoplanet : str
        The chosen target.
    archive : str, optional
        The exoplanet archive to pull data from. Options are:
        
        - 'eu'
        - 'nasa'
        The default is 'eu'.

    Returns
    -------
    None.

    '''
    temp = f'firefly/{exoplanet}'
    os.makedirs(temp, exist_ok=True)
    if archive == 'eu':
        if printing == False:
            with suppress_print():
                nan = _eu(exoplanet)
                nan = nan[7]
        elif printing == True:
            nan = _eu(exoplanet)
            nan = nan[7]
    elif archive == 'nasa':
        if printing == False:
            with suppress_print():
                nan = _nasa(exoplanet)
                nan = nan[7]
        elif printing == True:
            nan = _nasa(exoplanet)
            nan = nan[7]
    rmtree(temp)
    return nan