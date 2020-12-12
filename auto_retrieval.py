from transitfit import split_lightcurve_file, run_retrieval, calculate_logg
from lightkurve import search_lightcurvefile
from traceback import format_exc
from datetime import datetime, timedelta
from os import path, makedirs, remove, getcwd, walk, devnull
from smtplib import SMTP_SSL
from astropy.io import fits
from csv import DictWriter
from pandas import DataFrame, read_csv
from shutil import rmtree, move
from multiprocessing import Pool
import sys


class suppress_print():
    def __enter__(self):
        self.original_stdout = sys.stdout
        sys.stdout = open(devnull, 'w')
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self.original_stdout

def email(subject, body):
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

def auto_target(exoplanet):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    tess_filter_path = '/data/TESS_filter_path.csv'
    tess_filter = getcwd() + '/data/Filters/TESS_filter.csv'
    if not path.exists(tess_filter_path):
        cols = ['filter_idx', 'low_wl', 'high_wl']
        df = DataFrame(columns=cols)
        df = df.append([{'filter_idx': 0,
                         'low_wl': tess_filter}], ignore_index=True)
        df.to_csv(r'data/TESS_filter_path.csv', index=False, header=True)
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download EU Data
    makedirs('data', exist_ok=True)
    download_link = 'http://exoplanet.eu/catalog/csv'
    if not path.exists('data/eu.csv'):
        print('eu.csv does not exist, downloading.')
        df = read_csv(download_link)
        df.to_csv('data/eu.csv', index=False)
    ten_days_ago = datetime.now() - timedelta(days=10)
    filetime = datetime.fromtimestamp(path.getctime('data/eu.csv'))
    if filetime < ten_days_ago:
        print('eu.csv is 10 days old, updating.')
        df = read_csv(download_link)
        df.to_csv('data/eu.csv', index=False)
    else:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download MAST lightcurves
    lc = search_lightcurvefile(exoplanet, mission='TESS')
    print(lc)
    try:
        sector_list = lc .table .to_pandas()['sequence_number'] \
                         .drop_duplicates() .tolist()
    except Exception:
        sys.exit('Search result contains no data products for ' + exoplanet + '.')
    sector = sector_list
    for i in sector:
        sector = i
        sector_folder = 'Exoplanet/' + exoplanet + \
            '/TESS Sector ' + str(sector)
        if not path.exists(sector_folder):
            makedirs('Exoplanet/' + exoplanet + '/TESS Sector ' + str(sector),
                     exist_ok=True)
            lc = search_lightcurvefile(exoplanet, mission='TESS',
                                       sector=sector)
            print('\nDownloading MAST Lightcurve for TESS Sector ' +
                  str(sector) + '.')
            lc.download_all(download_dir='Exoplanet/' + exoplanet)
            source = 'Exoplanet/' + exoplanet + '/mastDownload/TESS/'
            if not path.exists(source):
                rmtree(sector_folder)
                sys.exit('No data products exist for the chosen TESS Sector.')
            files_in_dir = []
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Root, Directories, Files
            for r, d, f in walk(source):
                for item in f:
                    if '.fits' in item:
                        files_in_dir.append(path.join(r, item))
            destination = sector_folder + '/' + exoplanet + '.fits'
            for files in files_in_dir:
                if files.endswith(".fits"):
                    move(files, destination)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Extract Time series
            curvefile = 'Exoplanet/' + exoplanet + '/TESS Sector ' + \
                        str(sector) + '/split_curve_0.csv'
            if not path.exists(curvefile):
                fitsfile = 'Exoplanet/' + exoplanet + '/TESS Sector ' + \
                    str(sector) + '/' + exoplanet + '.fits'
                with fits.open(fitsfile) as TESS_fits:
                    time = TESS_fits[1].data['TIME'] + 2457000
                    time += TESS_fits[1].data['TIMECORR']
                    flux = TESS_fits[1].data['PDCSAP_FLUX']
                    flux_err = TESS_fits[1].data['PDCSAP_FLUX_ERR']
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                # Extract all light curves to a single csv file
                write_dict = []
                for i in range(len(time)):
                    write_dict.append({'Time': time[i], 'Flux': flux[i],
                                       'Flux err': flux_err[i]})
                csv_name = path.splitext(fitsfile)[0] + '.csv'

                with open(csv_name, 'w') as f:
                    columns = ['Time', 'Flux', 'Flux err']
                    writer = DictWriter(f, columns)
                    writer.writeheader()
                    writer.writerows(write_dict)
            else:
                pass
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Read in eu.csv
            col_subset_eu = [
                '# name',
                'orbital_period',
                'orbital_period_error_max',
                'semi_major_axis',
                'semi_major_axis_error_max',
                'radius',
                'radius_error_max',
                'eccentricity',
                'eccentricity_error_max',
                'inclination',
                'inclination_error_max',
                'tzero_tr',
                'tzero_tr_error_max',
                'star_teff',
                'star_teff_error_max',
                'star_radius',
                'star_radius_error_max',
                'star_mass',
                'star_mass_error_max',
                'star_metallicity',
                'star_metallicity_error_max']
            csv_file = read_csv('data/eu.csv', index_col='# name',
                                    usecols=col_subset_eu)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Pick Out Chosen Exoplanet Priors
            df = DataFrame(csv_file).loc[[exoplanet]]
            s = df.mean()
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
                     0.8 * radius_const * s.loc['radius'] / 
                                          s.loc['star_radius'],
                     1.2 * radius_const * s.loc['radius'] / 
                                          s.loc['star_radius'], 0]]
            repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                              'Input_A', 'Input_B', 'Filter'])
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Save the priors
            priors = 'Exoplanet/' + exoplanet + '/' + exoplanet + ' Priors.csv'
            if not path.exists(priors):
                repack.to_csv(priors, index=False, header=True)
            else:
                pass
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
                     0.8 * radius_const * s.loc['radius'] / 
                                          s.loc['star_radius'],
                     1.2 * radius_const * s.loc['radius'] / 
                                          s.loc['star_radius'], 0],
                    ['host_T', 'fixed', host_T[0], host_T[1], ''],
                    ['host_z', 'fixed', host_z[0], host_z[1], ''],
                    ['host_r', 'fixed', host_r[0], host_r[1], ''],
                    ['host_logg', 'fixed', host_logg[0], host_logg[1], '']]
            repack = DataFrame(cols, columns=['Parameter', 'Distribution',
                                              'Input_A', 'Input_B', 'Filter'])
            repack = repack.to_string(index=False)
            print('\nPriors generated from the EU Archive for ' + exoplanet +
                  ' are available to edit.\n')
            print(repack)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Split the Light curves
            if not path.exists(curvefile):
                t0, P = s.loc['tzero_tr'], s.loc['orbital_period']
                csvfile = 'Exoplanet/' + exoplanet + '/TESS Sector ' + \
                    str(sector) + '/' + exoplanet + '.csv'
                a = split_lightcurve_file(csvfile, t0=t0, P=P)
                print('\nA total of ' + str(len(a)) +
                      ' lightcurves were created.')
            else:
                pass
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Set the Data Paths
            cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
            df = DataFrame(columns=cols)
            # curves = len(a)
            curves = 1
            print()
            for i in range(curves):
                df = df.append([{'Path': getcwd() + '/Exoplanet/' + exoplanet +
                                 '/TESS Sector ' + str(sector) +
                                 '/split_curve_' + str(i) + '.csv'}],
                               ignore_index=True)
                df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
                df['Epochs'] = range(0, len(df))
            df.to_csv(r'data/data_paths.csv', index=False, header=True)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Cleanup
            try:
                rmtree('Exoplanet/' + exoplanet + '/mastDownload')
                remove('Exoplanet/' + exoplanet + '/' + '/TESS Sector ' +\
                      str(sector)
                       + '/' + exoplanet + '.csv')
                remove('Exoplanet/' + exoplanet + '/' + '/TESS Sector ' +\
                       str(sector)
                       + '/' + exoplanet + '.fits')
            except Exception:
                pass
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Paths to data, priors, and filter info:
            data = 'data/data_paths.csv'
            priors = 'Exoplanet/' + exoplanet + '/' + exoplanet + ' Priors.csv'
            filters = 'data/TESS_filter_path.csv'
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Output folders
            results_output_folder = 'Exoplanet/' + exoplanet + '/' + \
                '/TESS Sector ' + str(sector) + '/output_parameters'
            fitted_lightcurve_folder = 'Exoplanet/' + exoplanet + '/' + \
                '/TESS Sector ' + str(sector) + '/fitted_lightcurves'
            plot_folder = 'Exoplanet/' + exoplanet + '/' + '/TESS Sector ' +\
                str(sector) + '/plots'
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # Run the retrieval
            detrending = [['nth order', 2]]
            run_retrieval(data, priors, filters, detrending, host_T=host_T,
                          host_logg=host_logg, host_z=host_z, nlive = 1000,
                          host_r=host_r, dynesty_sample='rslice',
                          fitting_mode='folded', fit_ttv=True,
                          results_output_folder=results_output_folder,
                          final_lightcurve_folder=fitted_lightcurve_folder,
                          plot_folder=plot_folder)  

def email_handler(exoplanet_list):
    '''
    Automated version of target which inherits from auto_target. Sends an 
    email to transitfit.server@gmail.com upon an error or full completion of
    a target. Iteratively takes targets and employs TransitFit across each 
    sector for every exoplanet in the list given. Runs TransitFit for all 
    available split curves with the following options set:
        
        detrending = [['nth order', 2]]
        
        fitting_mode = 'folded'
        
        dynesty_sample = 'rslice'
        
        nlive = 1000
        
        fit_ttv = True

    Parameters
    ----------
    exoplanet_list : str
        A list of exoplanet targets.

    Returns
    -------
    A whole lot of data to science!

    '''
    for i in exoplanet_list:
        exoplanet_list = i
        try:
            # Printing suppressed within scope
            with suppress_print():
                auto_target(exoplanet_list)
            email('Success: ' + exoplanet_list, 
                  'Exoplanet: ' + exoplanet_list + '\n\n'
                  'A new target has been fully retrieved across ' +\
                  'all available TESS Sectors.')
        except KeyboardInterrupt:
            sys.exit('User terminated retrieval')
        except:
            trace_back = format_exc()
            email('Exception: ' + exoplanet_list, trace_back)
            pass    

def auto_retrieval(exoplanet_list, processes=4, chunksize=1):
    '''
    Automated version of target which inherits from auto_target. Sends an 
    email to transitfit.server@gmail.com upon an error or full completion of
    a target. Iteratively takes targets and employs TransitFit across each 
    sector for every exoplanet in the list given. Runs TransitFit for all 
    available split curves with the following options set:
        
        detrending = [['nth order', 2]]
        
        fitting_mode = 'folded'
        
        dynesty_sample = 'rslice'
        
        nlive = 1000
        
        fit_ttv = True

    Parameters
    ----------
    exoplanet_list : str
        A list of exoplanet targets.
    processes : int
        The number of processes to run in parallel. The default is 4.
    chunksize : int
        How many targets to assign to each process. The default is 1.
    
    Returns
    -------
    A whole lot of data to science!

    '''
    # Initialise Parallel Processing
    if __name__ == '__main__':
        with Pool(processes=processes) as pool:
            pool.map(email_handler, exoplanet_list, chunksize)


# Plans to read in from an external list file?
exoplanet_list = [
    ['WASP-91 b'], ['WASP-18 b'], ['WASP-43 b'], ['WASP-12 b'],
    ['WASP-126 b'], ['LHS 3844 b'], ['GJ 1252 b'], ['TOI-270 b']           
                 ]
auto_retrieval(exoplanet_list)
