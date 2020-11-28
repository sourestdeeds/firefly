import os
import pandas as pd
import numpy as np
from os import path
from datetime import datetime, timedelta


def target(exoplanet, curves = 1, dtype = 'eu'):
    '''
    Generates the priors and host star variables for a chosen target.
    Downloads exoplanet archives every 10 days and stores in /data.
    Target lightcurve files must be contained within a folder by the name of 
    its target.
    
    An example use with TransitFit is the following:
        
    Host Info
        exoplanet, curves, dtype = 'WASP-43 b', 1, 'eu'
        
        host_T, host_z, host_r, host_logg  = target(exoplanet, curves)
        
    Paths to data, priors, and filter info:
        data = 'data/data_paths.csv'
        
        priors = 'data/priors.csv'

    Parameters
    ----------
    exoplanet : string, example: 'WASP-43 b'
        The target exoplanet.
        
    curves : int, optional
        How many light curves to fit. Updates data paths for chosen target.
        Must be contained within the target folder, 
        ie 'WASP-43 b/split_curve_0.csv'.
        The default is 1.
        
    dtype : string, optional
        Allows for inputs 'nasa' or 'eu'. The default is 'eu'.
        
        EU : Data downloaded and stored in 'data/eu_data.csv'.
        http://exoplanet.eu/catalog/#
        
        NASA : Data downloaded and stored in 'data/nasa_data.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html

    Returns
    -------
    host_T : tuple or None
        The effective temperature of the host star, in Kelvin. 
    host_z : tuple or None
        The log_10 of the surface gravity of the host star, 
        with gravity measured in cm/s2. 
    host_r : tuple or None
        The metalicity of the host, given as a (value, uncertainty) pair.
    host_logg : tuple or None
        The host radius in Solar radii.
    data/data_paths.csv : file
        The locations of the light curves to fit.
    data/priors.csv : file
        The priors for the chosen target.
    data/eu.csv : file
        EU : Data downloaded and stored in 'data/eu_data.csv'.
        http://exoplanet.eu/catalog/#
    data/nasa.csv : file
        NASA : Data downloaded and stored in 'data/nasa_data.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download EU Data
    os.makedirs('data', exist_ok = True)
    
    if dtype == 'eu':
        download_link = 'http://exoplanet.eu/catalog/csv'
        
        if not os.path.exists('data/eu.csv'):
            print('eu.csv does not exist, downloading...')
            df = pd.read_csv(download_link)
            df.to_csv('data/eu.csv', index = False)
            
        ten_days_ago = datetime.now() - timedelta(days = 10)
        filetime = datetime.fromtimestamp(path.getctime('data/eu.csv'))
        if filetime < ten_days_ago:
            print('eu.csv is 10 days old, redownloading...')
            df = pd.read_csv(download_link)
            df.to_csv('data/eu.csv', index = False)
        else:
            print('eu.csv is recent.')
            pass
    # Download NASA Data
    else:
        download_link =                                                      \
            'https://exoplanetarchive.ipac.caltech.edu/'                    +\
            'TAP/sync?query=select+'                                        +\
            'pl_name,pl_orbper,pl_orbpererr1,pl_orbsmax,pl_orbsmaxerr1,'    +\
            'pl_radj,pl_radjerr1,pl_orbeccen,pl_orbeccenerr1,'              +\
            'st_teff,st_tefferr1,st_rad,st_raderr1,st_mass,st_masserr1,'    +\
            'st_met,st_meterr1,pl_tranmid,pl_tranmiderr1,'                  +\
            'pl_orbincl,pl_orbinclerr1,pl_orblper,pl_orblpererr1'           +\
            '+from+ps&format=csv'
        
        if not os.path.exists('data/nasa.csv'):
            print('nasa.csv does not exist, downloading...')
            df = pd.read_csv(download_link)
            df.to_csv('data/nasa.csv', index = False)
        ten_days_ago = datetime.now() - timedelta(days = 10)
        filetime = datetime.fromtimestamp(path.getctime('data/nasa.csv'))
        if filetime < ten_days_ago:
            print('nasa.csv is 10 days old, redownloading...')
            df = pd.read_csv(download_link)
            df.to_csv('data/nasa.csv', index = False)
        else:
            print('nasa.csv is recent.')
            pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = pd.DataFrame(columns = cols)  
    for i in range(curves):
        ############ UPDATE PATH TO LIGHTCURVES HERE #############
        # df = df.append([{'Path':'/Your/Path/Here/'
        df = df.append([{'Path':'/data/cmindoza/TransitFit/Planet/'
        ############ UPDATE PATH TO LIGHTCURVES HERE #############
                      +exoplanet+'/split_curve_'+str(i)+'.csv'}
                        ], ignore_index = True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    df.to_csv(r'data/data_paths.csv', index = False, header = True)
    print('Assigned data paths for '+exoplanet+'..')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Find Target
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
    if dtype == 'eu':
        csv_file = pd.read_csv('data/eu.csv', index_col = '# name',  
                               usecols = col_subset_eu)
    elif dtype == 'nasa':
        csv_file = pd.read_csv('data/nasa.csv', index_col = 'pl_name')  
                           #usecols = col_subset_nasa)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    df = pd.DataFrame(csv_file)
    df = df.loc[[exoplanet]] 
    s = df.mean() # Takes the mean of values if multiple entries
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    G = 6.674e-11
    m_sun = 1.989e30
    r_sun = 6.957e8
    if dtype == 'eu':
        # Host Logg Calc
        m = m_sun*s.loc['star_mass']
        err_m = m_sun*s.loc['star_mass_error_max']
        r = r_sun*s.loc['star_radius']
        err_r = r_sun*s.loc['star_radius_error_max']
        g = m * G * r ** -2 * 100
        err_g = G/(r ** 2) * np.sqrt(err_m ** 2 + 
                                     (4 * m**2 * err_r ** 2)/r**2) * 100
        logg = np.log10(g)
        err_logg = err_g / (g * np.log(10))
        
        host_T = (s.loc['star_teff'], s.loc['star_teff_error_max'])
        host_z = (s.loc['star_metallicity'], 
                      s.loc['star_metallicity_error_max'])
        host_r = (s.loc['star_radius'], s.loc['star_radius_error_max'])
        host_logg = ( logg, err_logg )
    elif dtype == 'nasa':
        # Host Logg Calc
        m = m_sun*s.loc['st_mass']
        err_m = m_sun*s.loc['st_masserr1']
        r = r_sun*s.loc['st_rad']
        err_r = r_sun*s.loc['st_raderr1']
        g = m * G * r ** -2 * 100
        err_g = G/(r ** 2) * np.sqrt(err_m ** 2 + 
                                     (4 * m**2 * err_r ** 2)/r**2) * 100
        logg = np.log10(g)
        err_logg = err_g / (g * np.log(10))
        
        host_T = (s.loc['st_teff'], s.loc['st_tefferr1'])
        host_z = (s.loc['st_met'], s.loc['st_meterr1'])
        host_r = (s.loc['st_rad'], s.loc['st_raderr1'])
        host_logg = ( logg, err_logg )
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Exoplanet Priors to TransitFit
    radius_const = 0.1027626851
    if dtype == 'eu':
        cols = [['P', 'gaussian', s.loc['orbital_period'], 
                     s.loc['orbital_period_error_max'], ''],
                ['t0', 'gaussian', s.loc['tzero_tr'], 
                     s.loc['tzero_tr_error_max'], ''],
                ['a', 'gaussian', s.loc['semi_major_axis'], 
                     s.loc['semi_major_axis_error_max'], ''],
                ['inc', 'gaussian', s.loc['inclination'], 
                     s.loc['inclination_error_max'], ''],
                ['rp', 'uniform',  
                     0.8*radius_const*s.loc['radius']/s.loc['star_radius'], 
                     1.2*radius_const*s.loc['radius']/s.loc['star_radius'], 0],
                ['ecc', 'fixed', s.loc['eccentricity'], 
                     s.loc['eccentricity_error_max'], ''],
                ['w', 'fixed', 90 , '' , '']]
    elif dtype == 'nasa':
        cols = [['P', 'gaussian', s.loc['pl_orbper'], 
                     s.loc['pl_orbpererr1'], ''],
                ['t0', 'gaussian', s.loc['pl_tranmid'], 
                     s.loc['pl_tranmiderr1'], ''],
                ['a', 'gaussian', s.loc['pl_orbsmax'], 
                     s.loc['pl_orbsmaxerr1'], ''],
                ['inc', 'gaussian', s.loc['pl_orbincl'], 
                     s.loc['pl_orbinclerr1'], ''],
                ['rp', 'uniform',  
                     0.8*radius_const*s.loc['pl_radj']/s.loc['st_rad'], 
                     1.2*radius_const*s.loc['pl_radj']/s.loc['st_rad'], 0],
                ['ecc', 'fixed', s.loc['pl_orbeccen'], 
                     s.loc['pl_orbeccenerr1'], ''],
                ['w', 'fixed', s.loc['pl_orblper'], 
                     s.loc['pl_orblpererr1'] , '']
               ]
    repack = pd.DataFrame(cols, columns = ['Parameter', 'Distribution', 
                                     'Input_A', 'Input_B', 'Filter'])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save the priors
    repack.to_csv(r'data/priors.csv', index = False, header = True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # For printing variables only
    if dtype == 'eu':
        cols = [['P', 'gaussian', s.loc['orbital_period'], 
                 s.loc['orbital_period_error_max'], ''],
                ['t0', 'gaussian', s.loc['tzero_tr'], 
                     s.loc['tzero_tr_error_max'], ''],
                ['a', 'gaussian', s.loc['semi_major_axis'], 
                     s.loc['semi_major_axis_error_max'], ''],
                ['inc', 'gaussian', s.loc['inclination'], 
                     s.loc['inclination_error_max'], ''],
                ['rp', 'uniform',  
                     0.8*radius_const*s.loc['radius']/s.loc['star_radius'], 
                     1.2*radius_const*s.loc['radius']/s.loc['star_radius'], 0],
                ['ecc', 'fixed', s.loc['eccentricity'], 
                     s.loc['eccentricity_error_max'], ''],
                    ['w', 'fixed', 90 , '' , ''],
                    ['host_T', 'fixed', host_T[0] , host_T[1] , ''],
                    ['host_z', 'fixed', host_z[0] , host_z[1] , ''],
                    ['host_r', 'fixed', host_r[0] , host_r[1] , ''],
                    ['host_logg', 'fixed', host_logg[0] , host_logg[1] , '']]
    elif dtype == 'nasa':
        cols = [['P', 'gaussian', s.loc['pl_orbper'], 
                 s.loc['pl_orbpererr1'], ''],
                ['t0', 'gaussian', s.loc['pl_tranmid'], 
                     s.loc['pl_tranmiderr1'], ''],
                ['a', 'gaussian', s.loc['pl_orbsmax'], 
                     s.loc['pl_orbsmaxerr1'], ''],
                ['inc', 'gaussian', s.loc['pl_orbincl'], 
                     s.loc['pl_orbinclerr1'], ''],
                ['rp', 'uniform',  
                     0.8*radius_const*s.loc['pl_radj']/s.loc['st_rad'], 
                     1.2*radius_const*s.loc['pl_radj']/s.loc['st_rad'], 0],
                ['ecc', 'fixed', s.loc['pl_orbeccen'], 
                     s.loc['pl_orbeccenerr1'], ''],
                ['w', 'fixed', s.loc['pl_orblper'], 
                     s.loc['pl_orblpererr1'] , ''],
                    ['host_T', 'fixed', host_T[0] , host_T[1] , ''],
                    ['host_z', 'fixed', host_z[0] , host_z[1] , ''],
                    ['host_r', 'fixed', host_r[0] , host_r[1] , ''],
                    ['host_logg', 'fixed', host_logg[0] , host_logg[1] , '']]
    repack = pd.DataFrame(cols, columns = ['Parameter', 'Distribution', 
                                     'Input_A', 'Input_B', 'Filter'])
    repack = repack.to_string(index = False)
    print('Assigned the following values for '+exoplanet+
                                                'and passing to TransitFit..')
    print(repack)
    return host_T, host_z, host_r, host_logg

# Host Info
exoplanet, curves, dtype = 'WASP-43 b', 27, 'eu'
host_T, host_z, host_r, host_logg  = target(exoplanet, curves, dtype)

# Paths to data, priors, and filter info:
data = 'data/data_paths.csv'
priors = 'data/priors.csv'
