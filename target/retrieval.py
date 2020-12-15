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
from functools import partial
import sys
import os



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
    except BaseException:
        # Continue on conn failure
        pass


def _TESS_filter():
    here = os.path.dirname(os.path.abspath(__file__))
    tess_filter_path = f'{here}/data/TESS_filter_path.csv'
    tess_filter = f'{here}/data/Filters/TESS_filter.csv'
    if not os.path.exists(tess_filter_path):
        cols = ['filter_idx', 'low_wl', 'high_wl']
        df = DataFrame(columns=cols)
        df = df.append([{'filter_idx': 0,
                         'low_wl': tess_filter}], ignore_index=True)
        df.to_csv(tess_filter_path, index=False, header=True)
    else:
        pass


def _eu(exoplanet):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download EU archive
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
    priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
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
    repack = repack.to_string(index=False)
    print(f'\nPriors generated from the EU Archive for {exoplanet}.\n')
    print(repack)
    return host_T, host_z, host_r, host_logg, t0, P


def nasa_full():
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
    df.to_csv('data/nasa_full.csv', index=False)
        
    
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
    priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
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
    repack = repack.to_string(index=False)
    print(f'\nPriors generated from the NASA Archive for {exoplanet}.\n')
    print(repack)
    return host_T, host_z, host_r, host_logg, t0, P, t14


def query(exoplanet, archive='eu'):
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
    try:
        sector_list = lc .table .to_pandas()['sequence_number'] \
                         .drop_duplicates() .tolist()
    except KeyError:
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    print(f'\n{lc}')
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


def retrieval(exoplanet, archive='eu', nlive=300, fit_ttv=False,
               detrending_list=[['nth order', 2]],
               dynesty_sample='rslice', fitting_mode='folded',
               limb_darkening_model='quadratic', ld_fit_method='independent',
               max_batch_parameters=25, batch_overlap=2, dlogz=None, 
               maxiter=None, maxcall=None, dynesty_bounding='multi', 
               normalise=True, detrend=True):
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

    nlive : int, optional
        The number of live points. The default is 1000.

    detrending :  optional
        Detrending. The default is [['nth order', 2]].

    dynesty_sample : str, optional
        Sampling method. The default is 'rslice'.

    fitting_mode : str, optional
        Fitting mode. The default is 'folded'.

    fit_ttv : boolean, optional
        Fit TTV. The default is False.

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
        df = df.append([{'Path': f'{os.getcwd()}/{sector_folder}' +
                         f'/split_curve_{str(i)}.csv'}],
                       ignore_index=True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    df.to_csv(r'data/data_paths.csv', index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = 'data/data_paths.csv'
    priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
    here = os.path.dirname(os.path.abspath(__file__))
    filters = f'{here}/data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'{sector_folder}/output_parameters'
    fitted_lightcurve_folder = f'{sector_folder}/fitted_lightcurves'
    plot_folder = f'{sector_folder}/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Run the retrieval
    run_retrieval(data, priors, filters, detrending_list=detrending_list,
                  host_T=host_T, host_logg=host_logg, host_z=host_z,
                  host_r=host_r, dynesty_sample=dynesty_sample,
                  fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                  results_output_folder=results_output_folder,
                  final_lightcurve_folder=fitted_lightcurve_folder,
                  plot_folder=plot_folder, nlive=nlive,
                  limb_darkening_model=limb_darkening_model, 
                  ld_fit_method=ld_fit_method,
                  max_batch_parameters=max_batch_parameters, 
                  batch_overlap=batch_overlap, dlogz=dlogz, 
                  maxiter=maxiter, maxcall=maxcall, 
                  dynesty_bounding=dynesty_bounding, 
                  normalise=normalise, detrend=detrend)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Cleanup
    try:
        rmtree(f'{sector_folder}/mastDownload')
        move(fitsfile, f'Exoplanet/{exoplanet}.fits')
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


def _retrieval(exoplanet, archive='eu', nlive=300, fit_ttv=False,
               detrending_list=[['nth order', 2]],
               dynesty_sample='rslice', fitting_mode='folded',
               limb_darkening_model='quadratic', ld_fit_method='independent',
               max_batch_parameters=25, batch_overlap=2, dlogz=None, 
               maxiter=None, maxcall=None, dynesty_bounding='multi', 
               normalise=True, detrend=True):
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
        print(f'\nA total of {str(curves)} lightcurves were created.\n')
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Set the Data Paths
        cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
        df = DataFrame(columns=cols)
        for i in range(curves):
            df = df.append([{'Path': f'{os.getcwd()}/{sector_folder}' +
                             f'/split_curve_{str(i)}.csv'}],
                           ignore_index=True)
            df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
            df['Epochs'] = range(0, len(df))
        df.to_csv(r'data/data_paths.csv', index=False, header=True)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Paths to data, priors, and filter info:
        data = 'data/data_paths.csv'
        priors = f'Exoplanet/{exoplanet}/{exoplanet} Priors.csv'
        here = os.path.dirname(os.path.abspath(__file__))
        filters = f'{here}/data/TESS_filter_path.csv'
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Output folders
        results_output_folder = f'{sector_folder}/output_parameters'
        fitted_lightcurve_folder = f'{sector_folder}/fitted_lightcurves'
        plot_folder = f'{sector_folder}/plots'
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Run the retrieval
        run_retrieval(data, priors, filters, detrending_list=detrending_list,
                      host_T=host_T, host_logg=host_logg, host_z=host_z,
                      host_r=host_r, dynesty_sample=dynesty_sample,
                      fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                      results_output_folder=results_output_folder,
                      final_lightcurve_folder=fitted_lightcurve_folder,
                      plot_folder=plot_folder, nlive=nlive,
                      limb_darkening_model=limb_darkening_model, 
                      ld_fit_method=ld_fit_method,
                      max_batch_parameters=max_batch_parameters, 
                      batch_overlap=batch_overlap, dlogz=dlogz, 
                      maxiter=maxiter, maxcall=maxcall, 
                      dynesty_bounding=dynesty_bounding, 
                      normalise=normalise, detrend=detrend)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Cleanup
        try:
            rmtree(f'{sector_folder}/mastDownload')
            move(fitsfile, f'Exoplanet/{exoplanet}.fits')
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


def _iterable_target(exoplanet_list, archive='eu', nlive=300,
                     detrending_list=[['nth order', 2]],
                     dynesty_sample='rslice', fitting_mode='folded', fit_ttv=False,
                     limb_darkening_model='quadratic', ld_fit_method='independent',
                     max_batch_parameters=25, batch_overlap=2, dlogz=None, 
                     maxiter=None, maxcall=None, dynesty_bounding='multi', 
                     normalise=True, detrend=True, email=False,
                     printing=False):
    for i, exoplanet in enumerate(exoplanet_list):
        try:
            # Printing suppressed within scope
            if printing == False:
                with suppress_print():
                    _retrieval(exoplanet, archive=archive, nlive=nlive,
                               detrending_list=detrending_list,
                               dynesty_sample=dynesty_sample,
                               fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                               limb_darkening_model=limb_darkening_model, 
                               ld_fit_method=ld_fit_method,
                               max_batch_parameters=max_batch_parameters, 
                               batch_overlap=batch_overlap, dlogz=dlogz, 
                               maxiter=maxiter, maxcall=maxcall, 
                               dynesty_bounding=dynesty_bounding, 
                               normalise=normalise, detrend=detrend)
            elif printing == True:
                _retrieval(exoplanet, archive=archive, nlive=nlive,
                           detrending_list=detrending_list,
                           dynesty_sample=dynesty_sample,
                           fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                           limb_darkening_model=limb_darkening_model, 
                           ld_fit_method=ld_fit_method,
                           max_batch_parameters=max_batch_parameters, 
                           batch_overlap=batch_overlap, dlogz=dlogz, 
                           maxiter=maxiter, maxcall=maxcall, 
                           dynesty_bounding=dynesty_bounding, 
                           normalise=normalise, detrend=detrend)
            if email == True:
                _email(f'Success: {exoplanet}',
                       f'Exoplanet: {exoplanet} \n\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.')
        except KeyboardInterrupt:
            sys.exit('User terminated retrieval')
        except BaseException:
            trace_back = format_exc()
            if email == True:
                _email(f'Exception: {exoplanet}', trace_back)
            else:
                print(trace_back)
            pass


def auto_retrieval(targets, processes=len(os.sched_getaffinity(0)) // 4,
                   archive='eu', nlive=300, detrending_list=[['nth order', 2]],
                   dynesty_sample='rslice', fitting_mode='folded', fit_ttv=False,
                   limb_darkening_model='quadratic', ld_fit_method='independent',
                   max_batch_parameters=25, batch_overlap=2, dlogz=None, 
                   maxiter=None, maxcall=None, dynesty_bounding='multi', 
                   normalise=True, detrend=True, email=False, printing=False):
    '''
    Automated version of retrieval. Sends an email to transitfit.server@gmail.com
    upon an error or full completion of a target. Iteratively takes targets and
    employs TransitFit across each TESS sector for every exoplanet in the list given.
    Runs TransitFit for all available split curves.
    
    Parameters
    ----------
    targets : str
        A list of exoplanet targets.
    processes : int, optional
        The number of processes to run in parallel. For UNIX, this default
        is the maximum available for the current process.
        The default is maximum available cores for the current process.
    archive: str, optional
        The exoplanet archive to use for priors. Supports 'eu' and 'nasa'.
        The default is 'eu'.
    nlive : int, optional
        The number of live points to use in the nested sampling retrieval.
        Default is 1000.
    detrending_list : array_like, shape (n_detrending_models, 2)
        A list of different detrending models. Each entry should consist
        of a method and a second parameter dependent on the method.
        Accepted methods are
            ['nth order', order]
            ['custom', function, [global fit indices, filter fit indices, epoch fit indices]]
            ['off', ]
        function here is a custom detrending function. TransitFit assumes
        that the first argument to this function is times and that all
        other arguments are single-valued - TransitFit cannot fit
        list/array variables. If 'off' is used, no detrending will be
        applied to the `LightCurve`s using this model.
        If a custom function is used, and some inputs to the function
        should not be fitted individually for each light curve, but should
        instead be shared either globally, within a given filter, or within
        a given epoch, the indices of where these fall within the arguments
        of the detrending function should be given as a list. If there are
        no indices to be given, then use an empty list: []
        e.g. if the detrending function is given by
            ```
            foo(times, a, b, c):
                # do something
            ```
        and a should be fitted globally, then the entry in the method_list
        would be ['custom', foo, [1], [], []].
    dynesty_sample : str, optional
        Method used to sample uniformly within the likelihood constraint,
        conditioned on the provided bounds. Unique methods available are:
        uniform sampling within the bounds('unif'), random walks with fixed
        proposals ('rwalk'), random walks with variable ("staggering")
        proposals ('rstagger'), multivariate slice sampling along preferred
        orientations ('slice'), "random" slice sampling along all
        orientations ('rslice'), "Hamiltonian" slices along random
        trajectories ('hslice'), and any callable function which follows
        the pattern of the sample methods defined in dynesty.sampling.
        'auto' selects the sampling method based on the dimensionality of
        the problem (from ndim). When ndim < 10, this defaults to 'unif'.
        When 10 <= ndim <= 20, this defaults to 'rwalk'. When ndim > 20,
        this defaults to 'hslice' if a gradient is provided and 'slice'
        otherwise. 'rstagger' and 'rslice' are provided as alternatives for
        'rwalk' and 'slice', respectively. Default is 'rslice'.
    fitting_mode : {`'auto'`, `'all'`, `'folded'`, `'batched'`}, optional
        The approach TransitFit takes towards limiting the number of parameters
        being simultaneously fitted. The available modes are:
        - `'auto'` : Will calculate the number of parameters required to fit
          all the data simulataneously. If this is less than max_parameters,
          will set to `'all'` mode, else will set to `'folded'` if at least one
          filter has at least 3 epochs in it. Otherwise will set to `'batched'`
        - `'all'` : Fits all parameters simultaneously, with no folding or
          batching of curves. Should be used with caution when fitting very
          large (~< 30) numbers of parameters.
        - `'folded'` : Useful for fitting curves with multiple epochs for each
          filter. TransitFit will fit each filter separately and produce a
          period-folded light curve for each filter, before fitting these
          simultaneously, using the `'batched'` approach if required.
        - `'batched'` : Useful for large numbers of light curves with
          relatively few shared filters, so `'folded'` loses large amounts of
          multi-epoch information. This mode splits the filters into sets of
          overlapping batches, runs each batch and uses the weighted means of
          each batch to produce a final result.
        Default is `'folded'`.
    fit_ttv : boolean, optional
        DESCRIPTION. The default is False.
    limb_darkening_model : str, optional
        The limb darkening model to use. Allowed models are
            - 'linear'
            - 'quadratic'
            - 'squareroot'
            - 'power2'
            - 'nonlinear'
        With the exception of the non-linear model, all models are constrained
        by the method in Kipping (2013), which can be found at
        https://arxiv.org/abs/1308.0009. Use `ldc_low_lim` and `ldc_high_lim`
        to control the behaviour of unconstrained coefficients.
        Default is 'quadratic'.
    ld_fit_method : {`'coupled'`, `'single'`, `'independent'`, `'off'`}, optional
        Determines the mode of fitting of limb darkening parameters. The
        available modes are:
            - `'coupled'` : all limb darkening parameters are fitted
              independently, but are coupled to a wavelength dependent
              model based on the host parameters through `ldkt`
            - `'single'` : LD parameters are still tied to a model, but
              only the first filter is actively fitted. The remaining
              filters are estimated based off the ratios given by ldtk for
              a host with the given parameters. This mode is useful for a
              large number of filters, as `'coupled'` or `'independent'`
              fitting will lead to much higher computation times.
            - `'independent'` : Each LD coefficient is fitted separately for
              each filter, with no coupling to the ldtk models.
            - `'off'` : Will use the fixed value provided in the input file
        Default is `'independent'`
    max_batch_parameters : int, optional
        The maximum number of parameters to use in a single retrieval.
        Default is 25.
    batch_overlap : int, optional
        The number of epochs to overlap in each batch. This will be adhered
        to where possible. Default is 2.
    dlogz : float, optional
        Retrieval iteration will stop when the estimated contribution of
        the remaining prior volume to the total evidence falls below this
        threshold. Explicitly, the stopping criterion is
        `ln(z + z_est) - ln(z) < dlogz`, where z is the current evidence
        from all saved samples and z_est is the estimated contribution from
        the remaining volume. The default is `1e-3 * (nlive - 1) + 0.01`.
    maxiter : int or `None`, optional
        The maximum number of iterations to run. If `None`, will
        continue until stopping criterion is reached. Default is `None`.
    maxcall : int or `None`, optional
        The maximum number of likelihood calls in retrieval. If None, will
        continue until stopping criterion is reached. Default is `None`.
    dynesty_bounding : {`'none'`, `'single'`, `'multi'`, `'balls'`, 'cubes'`}, optional
        The decomposition to use in sampling. Default is 'multi'
    normalise : bool, optional
        If True, will assume that the light curves have not been normalised and
        will fit normalisation constants within the retrieval. The range to
        fit normalisation constants c_n are automatically detected using
            ``1/f_min <= c_n <= 1/f_max``
        as the default range, where f_min and f_max are the minimum and maximum
        flux values for a given light curve. Default is True.
    detrend : bool, optional
        If True, will initialise detrending fitting. Default is True.
    email : bool, optional
        If True will send status emails. The default is False.
    printing : bool, optional
        If True will print outputs. The default is False.

    Returns
    -------
    A whole lot of data to science!

    '''
    exoplanet_list = []
    for i, exoplanet in enumerate(targets):
        exoplanet_list.append([exoplanet])
    func = partial(_iterable_target, 
                   archive=archive, nlive=nlive,
                   detrending_list=detrending_list, ld_fit_method=ld_fit_method,
                   dynesty_sample=dynesty_sample, fitting_mode=fitting_mode,
                   fit_ttv=fit_ttv, limb_darkening_model=limb_darkening_model, 
                   max_batch_parameters=max_batch_parameters, 
                   batch_overlap=batch_overlap, dlogz=dlogz, 
                   maxiter=maxiter, maxcall=maxcall, 
                   dynesty_bounding=dynesty_bounding, 
                   normalise=normalise, detrend=detrend, email=email,
                   printing=printing)
    with Pool(processes=processes) as pool:
        pool.map(func, exoplanet_list, chunksize=1)
