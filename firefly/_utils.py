from transitfit import split_lightcurve_file, run_retrieval, calculate_logg
from lightkurve import search_lightcurvefile
from traceback import format_exc
from datetime import datetime, timedelta
from smtplib import SMTP_SSL
from astropy.io import fits
from csv import DictWriter
from pandas import DataFrame, read_csv
from shutil import rmtree, move, make_archive
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
    os.makedirs('firefly/data', exist_ok=True)
    tess_filter_path = 'firefly/data/TESS_filter_path.csv'
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
    repack = repack.to_string(index=False)
    print(f'\nPriors generated from the EU Archive for {exoplanet}.\n')
    print(repack)
    return host_T, host_z, host_r, host_logg, t0, P, nan


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
    repack = repack.to_string(index=False)
    print(f'\nPriors generated from the NASA Archive for {exoplanet}.\n')
    print(repack)
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
                nan = nan[6]
        elif printing == True:
            nan = _eu(exoplanet)
            nan = nan[6]
    elif archive == 'nasa':
        if printing == False:
            with suppress_print():
                nan = _nasa(exoplanet)
                nan = nan[6]
        elif printing == True:
            nan = _nasa(exoplanet)
            nan = nan[6]
    rmtree(temp)
    return nan


def _fits(exoplanet, exo_folder):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Root, Directories, Files
    files_in_dir = []
    source = f'{exo_folder}/mastDownload/TESS/'
    for r, d, f in os.walk(source):
        for item in f:
            if '.fits' in item:
                files_in_dir.append(os.path.join(r, item))
    destination = f'{exo_folder}/{exoplanet}.fits'
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


def _retrieval(exoplanet, archive='eu', curve_sample=1, nlive=300, fit_ttv=False,
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
    except KeyError:
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    # Clear up previous sessions
    try:
        rmtree(f'firefly/{exoplanet}')
    except BaseException:
        pass
    curves_split, curves_delete = [], []
    for i, sector in enumerate(sector_list):
        exo_folder = f'firefly/{exoplanet}'
        os.makedirs(exo_folder, exist_ok=True)
        lc = search_lightcurvefile(exoplanet, mission='TESS',
                                   sector=sector)
        print(f'\nDownloading MAST Lightcurve for {exoplanet} -' +
              f' TESS Sector {str(sector)}.')
        lc.download_all(download_dir=exo_folder)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Extract all light curves to a single csv file
        fitsfile = _fits(exoplanet, exo_folder)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Download Archive
        with suppress_print():
            if archive == 'eu':
                host_T, host_z, host_r, host_logg, t0, P, nan = \
                                                    _eu(exoplanet)
            elif archive == 'nasa':
                host_T, host_z, host_r, host_logg, t0, P, t14, nan = \
                                                    _nasa(exoplanet)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Split the Light curves
        csvfile = f'{exo_folder}/{exoplanet}.csv'
        new_base_fname = f'sector_{sector}_split_curve'
        split_curves = split_lightcurve_file(csvfile, t0=t0, P=P, 
                                             new_base_fname=new_base_fname)
        curves = int(curve_sample * len(split_curves))
        if curves == 0:
            curves = 1
        curves_split.append(curves)
        curves_delete.append(len(split_curves))
        print(f'\nA total of {len(split_curves)} lightcurves were created.')
        print(f'\nA sample of {str(curves)} lightcurves from '
              f'TESS Sector {sector} will be used.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths
    data_path = f'{exo_folder}/data_paths.csv'
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns=cols)
    sector_curves = dict(zip(sector_list, curves_split))
    for sector, curves in sector_curves.items():
        for i in range(curves):
            df = df.append([{'Path': f'{os.getcwd()}/{exo_folder}' +
                             f'/sector_{sector}_split_curve_{i}.csv'}],
                           ignore_index=True)
            df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
            df['Epochs'] = range(0, len(df))
    print(f'\nIn total, {len(df)} lightcurves will be fitted across all'
          ' TESS Sectors.\n')
    df.to_csv(data_path, index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = data_path
    priors = f'{exo_folder}/{exoplanet} Priors.csv'
    filters = 'firefly/data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'{exo_folder}/output_parameters'
    fitted_lightcurve_folder = f'{exo_folder}/fitted_lightcurves'
    plot_folder = f'{exo_folder}/plots'
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
        rmtree(f'{exo_folder}/mastDownload')
        move(fitsfile, f'{exo_folder}.fits')
        os.remove(f'{exo_folder}.fits')
        os.remove(data)
        os.remove(csvfile)
        os.remove(priors)
        sector_curves = dict(zip(sector_list, curves_delete))
        for sector, curves in sector_curves.items():
            for i in range(curves):
                os.remove(f'{exo_folder}/sector_{sector}' +
                          f'_split_curve_{str(i)}.csv')
    except BaseException:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    make_archive(f'{exo_folder} {now}', format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{exoplanet}')
    rmtree(f'{exo_folder}')
    return f'{exo_folder} {now}.gz.tar'


def _iterable_target(exoplanet_list, archive='eu', curve_sample=1, nlive=300,
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
                    success = _retrieval(exoplanet, archive=archive, nlive=nlive,
                                         detrending_list=detrending_list,
                                         dynesty_sample=dynesty_sample,
                                         fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                                         limb_darkening_model=limb_darkening_model, 
                                         ld_fit_method=ld_fit_method,
                                         max_batch_parameters=max_batch_parameters, 
                                         batch_overlap=batch_overlap, dlogz=dlogz, 
                                         maxiter=maxiter, maxcall=maxcall, 
                                         dynesty_bounding=dynesty_bounding, 
                                         normalise=normalise, detrend=detrend,
                                         curve_sample=curve_sample)
            elif printing == True:
                success = _retrieval(exoplanet, archive=archive, nlive=nlive,
                                     detrending_list=detrending_list,
                                     dynesty_sample=dynesty_sample,
                                     fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                                     limb_darkening_model=limb_darkening_model, 
                                     ld_fit_method=ld_fit_method,
                                     max_batch_parameters=max_batch_parameters, 
                                     batch_overlap=batch_overlap, dlogz=dlogz, 
                                     maxiter=maxiter, maxcall=maxcall, 
                                     dynesty_bounding=dynesty_bounding, 
                                     normalise=normalise, detrend=detrend,
                                     curve_sample=curve_sample)
                print(f'\nData location: {success}\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.')
            if email == True:
                _email(f'Success: {exoplanet}',
                       f'Data location: {success} \n\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.')
        except KeyboardInterrupt:
            sys.exit('User terminated retrieval')
        except TypeError:
            trace_back = format_exc()
            if email == True:
                _email(f'Exception TypeError: {exoplanet}', trace_back)
            else:
                print(trace_back)
        except BaseException:
            trace_back = format_exc()
            if email == True:
                _email(f'Exception: {exoplanet}', trace_back)
            else:
                print(trace_back)
            pass
