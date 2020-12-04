from transitfit import split_lightcurve_file, run_retrieval
from lightkurve import search_lightcurvefile
from datetime import datetime, timedelta
from os import path, makedirs, remove, getcwd
from astropy.io import fits
from csv import DictWriter
from numpy import log, log10, sqrt
from pandas import DataFrame, read_csv
from shutil import rmtree
from sys import exit

def target(exoplanet, curves = 1, dtype = 'nasa'):
    '''
    A target data retriever for confirmed/candidate TESS exoplanets.
    Generates the priors and host star variables for a chosen target.
    Downloads exoplanet archives every 10 days and stores in /data.
    Target lightcurve files are downloaded from MAST, fits file is
    stored in Planet/exoplanet/exoplanet.fits.
    
    
    An example use with TransitFit is the following:
        
    Host Info
        exoplanet = 'WASP-43 b'
        
        host_T, host_z, host_r, host_logg  = target(exoplanet)
        
    Paths to data, priors, and filter info:
        data = 'data/data_paths.csv'
        
        priors = 'data/priors.csv'
        
    Outputs
        results_output_folder = 'Planet/'+exoplanet+'/output_parameters'
        
        fitted_lightcurve_folder = 'Planet/'+exoplanet+'/fitted_lightcurves'
    
        plot_folder = 'Planet/'+exoplanet+'/plots'

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
        Allows for inputs 'nasa' or 'eu'. The default is 'nasa'.
        
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
    Planet/exoplanet/exoplanet.fits : file
        MAST : Data downloaded and stored in 'Planet/exoplanet/exoplanet.fits'.
        https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download EU Data
    makedirs('data', exist_ok = True) 
    if dtype == 'eu':
        download_link = 'http://exoplanet.eu/catalog/csv'
        
        if not path.exists('data/eu.csv'):
            print('eu.csv does not exist, downloading...')
            df = read_csv(download_link)
            df.to_csv('data/eu.csv', index = False)
            
        ten_days_ago = datetime.now() - timedelta(days = 10)
        filetime = datetime.fromtimestamp(path.getctime('data/eu.csv'))
        if filetime < ten_days_ago:
            print('eu.csv is 10 days old, redownloading...')
            df = read_csv(download_link)
            df.to_csv('data/eu.csv', index = False)
        else:
            print('eu.csv is less than 10 days old.')
            pass
    # Download NASA Data
    else:
        download_link =                                                      \
            'https://exoplanetarchive.ipac.caltech.edu/'                    +\
            'TAP/sync?query=select+'                                        +\
            'pl_name,pl_orbper,pl_orbpererr1,pl_orbsmax,pl_orbsmaxerr1,'    +\
            'pl_radj,pl_radjerr1,pl_orbeccen,pl_orbeccenerr1,'              +\
            'st_teff,st_tefferr1,st_rad,st_raderr1,st_mass,st_masserr1,'    +\
            'st_met,st_meterr1,pl_tranmid,pl_tranmiderr1,pl_trandur,'       +\
            'pl_orbincl,pl_orbinclerr1,pl_orblper,pl_orblpererr1'           +\
            '+from+ps&format=csv'
        
        if not path.exists('data/nasa.csv'):
            print('nasa.csv does not exist, downloading...')
            df = read_csv(download_link)
            df.to_csv('data/nasa.csv', index = False)
        ten_days_ago = datetime.now() - timedelta(days = 10)
        filetime = datetime.fromtimestamp(path.getctime('data/nasa.csv'))
        if filetime < ten_days_ago:
            print('nasa.csv is 10 days old, redownloading...')
            df = read_csv(download_link)
            df.to_csv('data/nasa.csv', index = False)
        else:
            print('nasa.csv is less than 10 days old.')
            pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download MAST lightcurves
    try:
        if not path.exists('Planet/'+exoplanet+''):
            makedirs('Planet/'+exoplanet+'', exist_ok = True)   
            print('Downloading MAST Lightcurves...')
            lcfs = search_lightcurvefile(exoplanet, mission = 'TESS')         \
            .download_all(download_dir = 'Planet/'+exoplanet+'')              \
            .PDCSAP_FLUX.stitch().remove_nans()             
            lcfs.to_fits(path='Planet/'+exoplanet+'/'+exoplanet+'.fits', 
                                                               overwrite=True)
            rmtree('Planet/'+exoplanet+'/mastDownload')
        else:
            print('MAST Lightcurves already exist..')
            pass
    except Exception:
        rmtree('Planet/'+exoplanet+'')
        exit('Currently only TESS targets are supported.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    csvfile = 'Planet/'+exoplanet+'/'+exoplanet+'.csv'
    if not path.exists(csvfile):
        fitsfile = 'Planet/'+exoplanet+'/'+exoplanet+'.fits'
        print('Extracting the Time Series..')
        with fits.open(fitsfile) as TESS_fits:
            time = TESS_fits[1].data['TIME'] + 2457000
            flux = TESS_fits[1].data['FLUX']
            flux_err = TESS_fits[1].data['FLUX_ERR']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
        print('Writing csv file..')
        write_dict = []
        for i in range(len(time)):
            write_dict.append({'Time': time[i] , 'Flux':flux[i], 
                               'Flux err':flux_err[i]})   
        csv_name = path.splitext(fitsfile)[0] + '.csv'
    
        with open(csv_name, 'w') as f:
            columns = ['Time', 'Flux', 'Flux err']
            writer = DictWriter(f, columns)
            writer.writeheader()
            writer.writerows(write_dict)
        remove('Planet/'+exoplanet+'/'+exoplanet+'.fits')
    else:
        print('Lightcurves already extracted to csv..')
        pass
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
        csv_file = read_csv('data/eu.csv', index_col = '# name',  
                               usecols = col_subset_eu)
    elif dtype == 'nasa':
        csv_file = read_csv('data/nasa.csv', index_col = 'pl_name')  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    df = DataFrame(csv_file)
    df = df.loc[[exoplanet]] 
    s = df.mean() # Takes the mean of values if multiple entries
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign Host data to Transitfit
    G = 6.67430e-11
    m_sun = 1.98847e30
    r_sun = 6.957e8
    if dtype == 'eu':
        # Host Logg Calc
        m = m_sun*s.loc['star_mass']
        err_m = m_sun*s.loc['star_mass_error_max']
        r = r_sun*s.loc['star_radius']
        err_r = r_sun*s.loc['star_radius_error_max']
        g = m * G * r ** -2 * 100
        err_g = G/(r ** 2) * sqrt(err_m ** 2 + 
                                     (4 * m**2 * err_r ** 2)/r**2) * 100
        logg = log10(g)
        err_logg = err_g / (g * log(10))
        
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
        err_g = G/(r ** 2) * sqrt(err_m ** 2 + 
                                     (4 * m**2 * err_r ** 2)/r**2) * 100
        logg = log10(g)
        err_logg = err_g / (g * log(10))
        
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
                     1.2*radius_const*s.loc['radius']/s.loc['star_radius'], 0]]
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
                     1.2*radius_const*s.loc['pl_radj']/s.loc['st_rad'], 0]]
    repack = DataFrame(cols, columns = ['Parameter', 'Distribution', 
                                     'Input_A', 'Input_B', 'Filter'])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save the priors
    repack.to_csv(r'data/priors.csv', index = False, header = True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Split the Light curves
    curvefile = 'Planet/'+exoplanet+'/split_curve_0.csv'
    if not path.exists(curvefile):
        print('Splitting the light curves...')
        if dtype == 'eu':
            t0, P = s.loc['tzero_tr'] - 247000, s.loc['orbital_period']
            csvfile = 'Planet/'+exoplanet+'/'+exoplanet+'.csv'
            split_lightcurve_file(csvfile, t0, P)
        else:
            t0, P, t14 = s.loc['pl_tranmid'] - 247000, s.loc['pl_orbper'], \
                                s.loc['pl_trandur']*24
            csvfile = 'Planet/'+exoplanet+'/'+exoplanet+'.csv'
            split_lightcurve_file(csvfile, t0, P, t14)
    else:
        print('Lightcurves already split..')
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns = cols)  
    for i in range(curves):
        df = df.append([{'Path':getcwd()
                         +'/Planet/'+exoplanet+'/split_curve_'+str(i)+'.csv'}
                        ], ignore_index = True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    df.to_csv(r'data/data_paths.csv', index = False, header = True)
    print('Assigned data paths for '+exoplanet+'..')
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
                ['host_T', 'fixed', host_T[0] , host_T[1] , ''],
                ['host_z', 'fixed', host_z[0] , host_z[1] , ''],
                ['host_r', 'fixed', host_r[0] , host_r[1] , ''],
                ['host_logg', 'fixed', host_logg[0] , host_logg[1] , '']]
    repack = DataFrame(cols, columns = ['Parameter', 'Distribution', 
                                     'Input_A', 'Input_B', 'Filter'])
    repack = repack.to_string(index = False)
    print('Assigned the following values for '+exoplanet+
                                                ' and passing to TransitFit..')
    print(repack)
    return host_T, host_z, host_r, host_logg

def main(exoplanet, curves):
    host_T, host_z, host_r, host_logg = target(exoplanet, curves)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Override if data not present - Guess!
    # host_T = ( , )
    # host_z = (0.1 , 0.1)
    # host_r = ( , )
    # host_logg = ( , )
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = 'data/data_paths.csv'
    priors = 'data/priors.csv'
    filters = 'data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    detrending = [['nth order', 2]]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Outputs
    results_output_folder = 'Planet/'+exoplanet+'/output_parameters'
    fitted_lightcurve_folder = 'Planet/'+exoplanet+'/fitted_lightcurves'
    plot_folder = 'Planet/'+exoplanet+'/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Run the retrieval
    run_retrieval(data, priors, filters, detrending, host_T=host_T,
                  host_logg = host_logg, host_z = host_z, 
                  host_r=host_r[0], dynesty_sample = 'rslice',
                  fitting_mode = 'folded', fit_ttv = True,
                  results_output_folder = results_output_folder,
                  final_lightcurve_folder = fitted_lightcurve_folder,
                  plot_folder = plot_folder)

main('WASP-43 b', 1)
