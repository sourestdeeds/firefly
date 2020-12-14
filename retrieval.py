from transitfit import split_lightcurve_file, run_retrieval, calculate_logg
from lightkurve import search_lightcurvefile
from traceback import format_exc
from datetime import datetime, timedelta
from smtplib import SMTP_SSL
from astropy.io import fits
from csv import DictWriter
from pandas import DataFrame, read_csv
from shutil import rmtree, move, make_archive
from multiprocessing import Pool, cpu_count
import sys, os


class suppress_print():
    def __enter__(self):
        self.original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self.original_stdout

def _email(subject, body):
    username = 'transitfit.server@gmail.com'
    password = 'yvoq efzi dcib dmbm'
    sent_from = username
    to = ['transitfit.server@gmail.com']
    message = f'Subject: {subject}\n\n{body}'
    try:
        server = SMTP_SSL('smtp.gmail.com', 465)
        server.ehlo()
        server.login(username, password)
        server.sendmail(sent_from, to, message)
        server.close()
    except:
        # Continue on conn failure
        pass

def _TESS_filter():
    tess_filter_path = '/data/TESS_filter_path.csv'
    tess_filter = f'{os.getcwd()}/data/Filters/TESS_filter.csv'
    if not os.path.exists(tess_filter_path):
        cols = ['filter_idx', 'low_wl', 'high_wl']
        df = DataFrame(columns=cols)
        df = df.append([{'filter_idx': 0,
                         'low_wl': tess_filter}], ignore_index=True)
        df.to_csv(r'data/TESS_filter_path.csv', index=False, header=True)
    else:
        pass

def _eu(exoplanet):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download EU exo archive
    os.makedirs('data', exist_ok=True)
    download_link = 'http://exoplanet.eu/catalog/csv'
    if not os.path.exists('data/eu.csv'):
        print('eu.csv does not exist, downloading.')
        df = read_csv(download_link)
        df.to_csv('data/eu.csv', index=False)
    ten_days_ago = datetime.now() - timedelta(days=10)
    filetime = datetime.fromtimestamp(os.path.getctime('data/eu.csv'))
    if filetime < ten_days_ago:
        print('eu.csv is 10 days old, updating.')
        df = read_csv(download_link)
        df.to_csv('data/eu.csv', index=False)
    else:
        pass
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
    exo_archive = read_csv('data/eu.csv', index_col='# name',
                            usecols=col_subset_eu)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    df = DataFrame(exo_archive).loc[[exoplanet]]
    s = df.mean()
    t0, P = s.loc['tzero_tr'], s.loc['orbital_period']
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
             s.loc['orbital_period_error_max'], ''],
            ['t0', 'gaussian', s.loc['tzero_tr'],
             s.loc['tzero_tr_error_max'], ''],
            ['a', 'gaussian', s.loc['semi_major_axis'],
             s.loc['semi_major_axis_error_max'], ''],
            ['inc', 'gaussian', s.loc['inclination'],
             s.loc['inclination_error_max'], ''],
            ['rp', 'uniform',
             0.8 * radius_const * s.loc['radius'] / s.loc['star_radius'],
             1.2 * radius_const * s.loc['radius'] / s.loc['star_radius'], 0]]
    priors_csv = DataFrame(cols, columns=['Parameter', 'Distribution',
                          'Input_A', 'Input_B', 'Filter'])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save the priors
    priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
    if not os.path.exists(priors):
        priors_csv.to_csv(priors, index=False, header=True)
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # For printing variables only
    cols = [['P', 'gaussian', s.loc['orbital_period'],
                 s.loc['orbital_period_error_max'], ''],
                ['t0', 'gaussian', s.loc['tzero_tr'],
                 s.loc['tzero_tr_error_max'], ''],
                ['a', 'gaussian', s.loc['semi_major_axis'],
                 s.loc['semi_major_axis_error_max'], ''],
                ['inc', 'gaussian', s.loc['inclination'],
                 s.loc['inclination_error_max'], ''],
                ['rp', 'uniform',
                 0.8 * radius_const * s.loc['radius'] / s.loc['star_radius'],
                 1.2 * radius_const * s.loc['radius'] / s.loc['star_radius'], 0],
                ['host_T', 'fixed', host_T[0], host_T[1], ''],
                ['host_z', 'fixed', host_z[0], host_z[1], ''],
                ['host_r', 'fixed', host_r[0], host_r[1], ''],
                ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input_A', 'Input_B', 'Filter'])
    repack = repack.to_string(index=False)
    print(f'\nPriors generated from the EU Archive for {exoplanet}.\n')
    print(repack)
    return host_T, host_z, host_r, host_logg, t0, P

def _nasa(exoplanet):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download NASA archive
    download_link =  \
        'https://exoplanetarchive.ipac.caltech.edu/' +\
        'TAP/sync?query=select+' +\
        'pl_name,pl_orbper,pl_orbpererr1,pl_orbsmax,pl_orbsmaxerr1,' +\
        'pl_radj,pl_radjerr1,pl_orbeccen,pl_orbeccenerr1,disc_facility,' +\
        'st_teff,st_tefferr1,st_rad,st_raderr1,st_mass,st_masserr1,' +\
        'st_met,st_meterr1,pl_tranmid,pl_tranmiderr1,pl_trandur,' +\
        'pl_orbincl,pl_orbinclerr1,pl_orblper,pl_orblpererr1,rowupdate' +\
        '+from+ps&format=csv'
    if not os.path.exists('data/nasa.csv'):
        print('nasa.csv does not exist, downloading.')
        df = read_csv(download_link)
        df.to_csv('data/nasa.csv', index=False)
    ten_days_ago = datetime.now() - timedelta(days=10)
    filetime = datetime.fromtimestamp(os.path.getctime('data/nasa.csv'))
    if filetime < ten_days_ago:
        print('nasa.csv is 10 days old, updating.')
        df = read_csv(download_link)
        df.to_csv('data/nasa.csv', index=False)
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Read in nasa.csv
    exo_archive = read_csv('data/nasa.csv', index_col='pl_name')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pick Out Chosen Exoplanet Priors
    df = DataFrame(exo_archive).loc[[exoplanet]]
    s = df.mean()
    t0, P, t14 = s.loc['pl_tranmid'], s.loc['pl_orbper'], \
            s.loc['pl_trandur'] * 24 * 60
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
             s.loc['pl_orbpererr1']*1e5, ''],
            ['t0', 'gaussian', s.loc['pl_tranmid'],
             s.loc['pl_tranmiderr1']*100, ''],
            ['a', 'gaussian', s.loc['pl_orbsmax'],
             s.loc['pl_orbsmaxerr1'], ''],
            ['inc', 'gaussian', s.loc['pl_orbincl'],
             s.loc['pl_orbinclerr1'], ''],
            ['rp', 'uniform',
             0.8 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'],
             1.2 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'], 0]]
    priors_csv = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input_A', 'Input_B', 'Filter'])
    priors_csv = DataFrame(cols, columns=['Parameter', 'Distribution',
                          'Input_A', 'Input_B', 'Filter'])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save the priors
    priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
    if not os.path.exists(priors):
        priors_csv.to_csv(priors, index=False, header=True)
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # For printing variables only
    cols = [['P', 'gaussian', s.loc['pl_orbper'],
                 s.loc['pl_orbpererr1'], ''],
                ['t0', 'gaussian', s.loc['pl_tranmid'],
                 s.loc['pl_tranmiderr1'], ''],
                ['a', 'gaussian', s.loc['pl_orbsmax'],
                 s.loc['pl_orbsmaxerr1'], ''],
                ['inc', 'gaussian', s.loc['pl_orbincl'],
                 s.loc['pl_orbinclerr1'], ''],
                ['rp', 'uniform',
                 0.8 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'],
                 1.2 * radius_const * s.loc['pl_radj'] / s.loc['st_rad'], 0],
                ['host_T', 'fixed', host_T[0], host_T[1], ''],
                ['host_z', 'fixed', host_z[0], host_z[1], ''],
                ['host_r', 'fixed', host_r[0], host_r[1], ''],
                ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]
    repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                      'Input_A', 'Input_B', 'Filter'])
    repack = repack.to_string(index=False)
    print(f'\nPriors generated from the NASA Archive for {exoplanet}.\n')
    print(repack)
    return host_T, host_z, host_r, host_logg, t0, P, t14

def query(exoplanet, archive = 'eu'):
    '''
    Performs a query for prior information and data products from MAST

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
    temp = f'Exoplanet/{exoplanet}'
    os.makedirs(temp, exist_ok=True)
    lc = search_lightcurvefile(exoplanet, mission='TESS')
    print(lc)
    if archive == 'eu':
        _eu(exoplanet)
    elif archive == 'nasa':
        _nasa(exoplanet)
    rmtree(temp)

def _fits(exoplanet, sector_folder):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Root, Directories, Files
    files_in_dir = []
    source = f'{sector_folder}/mastDownload/TESS/'
    for r, d, f in os.walk(source):
        for item in f:
            if '.fits' in item:
                files_in_dir.append(os.path.join(r, item))
    destination = f'{sector_folder}/{exoplanet}.fits'
    for files in files_in_dir:
        if files.endswith(".fits"):
            move(files, destination)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    fitsfile = destination
    with fits.open(fitsfile) as TESS_fits:
        time = TESS_fits[1].data['TIME'] + 2457000
        time += TESS_fits[1].data['TIMECORR']
        flux = TESS_fits[1].data['PDCSAP_FLUX']
        flux_err = TESS_fits[1].data['PDCSAP_FLUX_ERR']
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
    write_dict = []
    for i in range(len(time)):
        write_dict.append({'Time': time[i], 'Flux': flux[i],
                           'Flux err': flux_err[i]})
    csv_name = os.path.splitext(fitsfile)[0] + '.csv'
    with open(csv_name, 'w') as f:
        columns = ['Time', 'Flux', 'Flux err']
        writer = DictWriter(f, columns)
        writer.writeheader()
        writer.writerows(write_dict)
    return fitsfile

def retrieval(exoplanet, archive='eu'):
    '''
    A target data retriever for confirmed/candidate TESS exoplanets.
    Generates the priors and host star variables for a chosen target.
    Downloads exoplanet archives every 10 days and stores in /data.
    Target lightcurve files are downloaded from MAST, fits file is
    stored in Planet/exoplanet/exoplanet.fits.


    An example use with TransitFit is the following:

    Host Info
        exoplanet = 'WASP-43 b'

        retrieval(exoplanet)

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

    archive : string, optional
        Allows for inputs 'nasa' or 'eu'. The default is 'nasa'.

        EU : Data downloaded and stored in 'data/eu_data.csv'.
        http://exoplanet.eu/catalog/#

        NASA : Data downloaded and stored in 'data/nasa_data.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html

    Returns
    -------
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
    # Check inputs are sensible
    if not archive == 'eu' or archive == 'nasa':
        sys.exit('Archive data options for dtype are: \'eu\' or \'nasa\'')
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    _TESS_filter()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download MAST lightcurves
    lc = search_lightcurvefile(exoplanet, mission='TESS')
    print(lc)
    try:
        sector_list = lc .table .to_pandas()['sequence_number'] \
                         .drop_duplicates() .tolist()
    except Exception:
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    sector = int(input('Enter which TESS Sector you' +
                       ' would like to download: '))
    sector_folder = f'Exoplanet/{exoplanet}/TESS Sector {str(sector)}'
    os.makedirs(sector_folder, exist_ok=True)
    lc = search_lightcurvefile(exoplanet, mission='TESS',
                                 sector=sector)
    print(f'\nDownloading MAST Lightcurve for TESS Sector {str(sector)}.')
    lc.download_all(download_dir=sector_folder)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
    fitsfile = _fits(exoplanet, sector_folder)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download Archive
    if archive == 'eu':
        host_T, host_z, host_r, host_logg, t0, P = _eu(exoplanet)
    elif archive == 'nasa':
        host_T, host_z, host_r, host_logg, t0, P, t14 = _nasa(exoplanet)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Split the Light curves
    csvfile = f'{sector_folder}/{exoplanet}.csv'
    split_curves = split_lightcurve_file(csvfile, t0=t0, P=P)
    print(f'\nA total of {str(len(split_curves))} lightcurves were created.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns=cols)
    curves = int(input('Enter how many lightcurves you wish to fit: '))
    print()
    for i in range(curves):
        df = df.append([{'Path': f'{os.getcwd()}/{sector_folder}' +\
                         f'/split_curve_{str(i)}.csv'}],
                       ignore_index=True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    df.to_csv(r'data/data_paths.csv', index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = 'data/data_paths.csv'
    priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
    filters = 'data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'Exoplanet/{exoplanet}' +\
        f'/TESS Sector {str(sector)}/output_parameters'
    fitted_lightcurve_folder = f'Exoplanet/{exoplanet}' +\
        f'/TESS Sector {str(sector)}/fitted_lightcurves'
    plot_folder = f'Exoplanet/{exoplanet}/TESS Sector {str(sector)}/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Run the retrieval
    detrending = [['nth order', 2]]
    run_retrieval(data, priors, filters, detrending, host_T=host_T,
                  host_logg=host_logg, host_z=host_z, nlive = 1000,
                  host_r=host_r, dynesty_sample='rslice',
                  fitting_mode='folded', fit_ttv=False,
                  results_output_folder=results_output_folder,
                  final_lightcurve_folder=fitted_lightcurve_folder,
                  plot_folder=plot_folder)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Cleanup
    try:
        rmtree(f'{sector_folder}/mastDownload')
        move(fitsfile,f'Exoplanet/{exoplanet}.fits')
        os.remove(f'Exoplanet/{exoplanet}.fits')
        os.remove(csvfile)
        os.remove(priors)
        for i in range(len(split_curves)):
            os.remove(f'{sector_folder}/split_curve_{str(i)}.csv')
    except Exception:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    make_archive(f'Exoplanet/{exoplanet} {now}', format='gztar', 
                 root_dir=f'{os.getcwd()}/Exoplanet/',
                 base_dir=f'{exoplanet}')
    rmtree(f'Exoplanet/{exoplanet}')

def _retrieval(exoplanet, archive='eu'):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    _TESS_filter()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download MAST lightcurves
    lc = search_lightcurvefile(exoplanet, mission='TESS')
    print(lc)
    try:
        sector_list = lc .table .to_pandas()['sequence_number'] \
                         .drop_duplicates() .tolist()
    except Exception:
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    # Clear up previous sessions
    try:
        rmtree(f'Exoplanet/{exoplanet}')
    except Exception:
        pass
    for i, sector in enumerate(sector_list):
        sector_folder = f'Exoplanet/{exoplanet}/TESS Sector {str(sector)}'
        os.makedirs(sector_folder, exist_ok=True)
        lc = search_lightcurvefile(exoplanet, mission='TESS',
                                   sector=sector)
        print(f'\nDownloading MAST Lightcurve for TESS Sector {str(sector)}.')
        lc.download_all(download_dir=sector_folder)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Extract all light curves to a single csv file
        fitsfile = _fits(exoplanet, sector_folder)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Download Archive
        if archive == 'eu':
            host_T, host_z, host_r, host_logg, t0, P = _eu(exoplanet)
        elif archive == 'nasa':
            host_T, host_z, host_r, host_logg, t0, P, t14 = _nasa(exoplanet)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Split the Light curves
        csvfile = f'{sector_folder}/{exoplanet}.csv'
        split_curves = split_lightcurve_file(csvfile, t0=t0, P=P)
        curves = len(split_curves)
        print(f'\nA total of {str(curves)} lightcurves were created.')
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Set the Data Paths
        cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
        df = DataFrame(columns=cols)
        for i in range(curves):
            df = df.append([{'Path': f'{os.getcwd()}/{sector_folder}' +\
                             f'/split_curve_{str(i)}.csv'}],
                           ignore_index=True)
            df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
            df['Epochs'] = range(0, len(df))
        df.to_csv(r'data/data_paths.csv', index=False, header=True)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Paths to data, priors, and filter info:
        data = 'data/data_paths.csv'
        priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
        filters = 'data/TESS_filter_path.csv'
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Output folders
        results_output_folder = f'Exoplanet/{exoplanet}' +\
            f'/TESS Sector {str(sector)}/output_parameters'
        fitted_lightcurve_folder = f'Exoplanet/{exoplanet}' +\
            f'/TESS Sector {str(sector)}/fitted_lightcurves'
        plot_folder = f'Exoplanet/{exoplanet}/TESS Sector {str(sector)}/plots'
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Run the retrieval
        detrending = [['nth order', 2]]
        run_retrieval(data, priors, filters, detrending, host_T=host_T,
                      host_logg=host_logg, host_z=host_z, nlive = 1000,
                      host_r=host_r, dynesty_sample='rslice',
                      fitting_mode='folded', fit_ttv=False,
                      results_output_folder=results_output_folder,
                      final_lightcurve_folder=fitted_lightcurve_folder,
                      plot_folder=plot_folder)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Cleanup
        try:
            rmtree(f'{sector_folder}/mastDownload')
            move(fitsfile,f'Exoplanet/{exoplanet}.fits')
            os.remove(f'Exoplanet/{exoplanet}.fits')
            os.remove(csvfile)
            os.remove(priors)
            for i in range(len(split_curves)):
                os.remove(f'{sector_folder}/split_curve_{str(i)}.csv')
        except Exception:
            pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    make_archive(f'Exoplanet/{exoplanet} {now}', format='gztar', 
                 root_dir=f'{os.getcwd()}/Exoplanet/',
                 base_dir=f'{exoplanet}')
    rmtree(f'Exoplanet/{exoplanet}')
 
def _iterable_target_eu(exoplanet_list):
    for i, exoplanet in enumerate(exoplanet_list):
        try:
            # Printing suppressed within scope
            with suppress_print():
                _retrieval(exoplanet, archive='eu')
            _email(f'Success: {exoplanet}', 
                  f'Exoplanet: {exoplanet} \n\n'
                  'A new target has been fully retrieved across ' +\
                  'all available TESS Sectors.')
        except KeyboardInterrupt:
            sys.exit('User terminated retrieval')
        except:
            trace_back = format_exc()
            _email(f'Exception: {exoplanet}', trace_back)
            pass

def _iterable_target_nasa(exoplanet_list):
    for i, exoplanet in enumerate(exoplanet_list):
        try:
            # Printing suppressed within scope
            with suppress_print():
                _retrieval(exoplanet, archive='nasa')
            _email(f'Success: {exoplanet}', 
                  f'Exoplanet: {exoplanet} \n\n'
                  'A new target has been fully retrieved across ' +\
                  'all available TESS Sectors.')
        except KeyboardInterrupt:
            sys.exit('User terminated retrieval')
        except:
            trace_back = format_exc()
            _email(f'Exception: {exoplanet}', trace_back)
            pass    

def auto_retrieval(exoplanet_list, processes=len(os.sched_getaffinity(0)),
                   archive='eu'):
    '''
    Automated version of retrieval. Sends an email to transitfit.server@gmail.com 
    upon an error or full completion of a target. Iteratively takes targets and 
    employs TransitFit across each TESS sector for every exoplanet in the list given. 
    Runs TransitFit for all available split curves with the following options set:
        
        detrending = [['nth order', 2]]
        
        fitting_mode = 'folded'
        
        dynesty_sample = 'rslice'
        
        nlive = 1000
        
        fit_ttv = False

    Parameters
    ----------
    exoplanet_list : str
        A list of exoplanet targets.
    archive: str, optional
        The exoplanet archive to use for priors. Supports 'eu' and 'nasa'.
        The default is 'eu'.
    processes : int, optional
        The number of processes to run in parallel. For UNIX, this default 
        is the maximum available for the current process.
        The default is maximum available cores for the current process.
    chunksize : int, optional
        How many targets to assign to each process. The default is 1.
    
    Returns
    -------
    A whole lot of data to science!

    '''
    with Pool(processes=processes) as pool:
        if archive == 'eu':
            pool.map(_iterable_target_eu, exoplanet_list, chunksize=1)
        elif archive == 'nasa':
            pool.map(_iterable_target_nasa, exoplanet_list, chunksize=1)

# Example use
'''
query('WASP_43 b')

retrieval('WASP-43 b')

exoplanet_list = (['WASP-43 b'], ['WASP-18 b'], ['WASP-91 b'], ['WASP-12 b'],
                  ['WASP-126 b'], ['LHS 3844 b'], ['GJ 1252 b'], ['TOI-270 b'])

auto_retrieval(exoplanet_list)

'''
    
