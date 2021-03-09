#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The backend for auto_retrieval.

@author: Steven Charles-Mindoza
"""

from ._archive import priors, _tic, _load_csv, _search
from ._plot import oc_fold, density_scatter

from transitfit import split_lightcurve_file, run_retrieval
from astroquery.mast import Observations as obs
from datetime import datetime
from tabulate import tabulate
from astropy.table import Table
from pandas import DataFrame, read_csv, Categorical
from shutil import rmtree, make_archive
from natsort import natsorted
from math import ceil
import numpy as np
import sys
import os
import random



def _TESS_filter():
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs('firefly/data', exist_ok=True)
    tess_filter_path = 'firefly/data/TESS_filter_path.csv'
    try:
        os.remove(tess_filter_path)
    except:
        pass
    tess_filter = f'{here}/data/Filters/TESS_filter.csv'
    cols = ['filter_idx', 'low_wl', 'high_wl']
    df = DataFrame(columns=cols)
    df = df.append([{'filter_idx': 0,
                     'low_wl': tess_filter}], ignore_index=True)
    df.to_csv(tess_filter_path, index=False, header=True)


def _alias(exoplanet):
    alias = read_csv('https://exoplanetarchive.ipac.caltech.edu/cgi-bin/' +\
        f'nstedAPI/nph-nstedAPI?table=aliastable&objname={exoplanet[:-2]}')
    check_wasp = alias['aliasdis'][1]
    if '1SWASP' in check_wasp:
        return check_wasp
    else:
        return None

def mast(exoplanet):
    _load_csv()
    highest, ratios = _search(exoplanet)
    exoplanet = highest[0]
    tic_id = _tic(exoplanet).replace('TIC ', '')
    print(f'\nSearching MAST for {exoplanet} (TIC {tic_id}).')
    search = obs.query_criteria(dataproduct_type=['timeseries'],
                                project='TESS',
                                target_name=tic_id
                                ).to_pandas()
    data = search[['obs_id', 'target_name', 't_exptime',
               'provenance_name', 'project']]
    data = data.rename(columns={"obs_id": "Product",
                                "target_name": "TIC ID",
                                "t_exptime": "Cadence",
                                "provenance_name": "HLSP",
                                "project": "Mission"})
    data['Product'] = \
                Categorical(data['Product'],
                ordered=True,
                categories=natsorted(data['Product'].unique()))
    data = data.sort_values('Product')
    print(f'\nQuery from MAST returned {len(search)} '
          f'data products for {exoplanet} (TIC {tic_id}).\n')
    return print(tabulate(data, tablefmt='psql', showindex=False, headers='keys'))


def discover():
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs('firefly/discover', exist_ok=True)
    upper_ecliptic = f'{here}/data/Candidates/19_sector_candidates.csv'
    #upper_ecliptic = '19_sector_candidates.csv'
    tic_ids = read_csv(upper_ecliptic)['TIC ID'].tolist()
    for i, tic_id in enumerate(tic_ids):
        _discover_search(tic_id)
        
        
def _discover_search(tic_id):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # MAST Search
    search = obs.query_criteria(dataproduct_type=['timeseries'],
                                project='TESS',
                                provenance_name='SPOC',
                                t_exptime=120,
                                target_name=tic_id
                                ).to_pandas()
    data = search[['obs_id', 'target_name', 'dataURL', 't_exptime',
                   'provenance_name', 'project']]
    data['dataURL'] = ['https://mast.stsci.edu/api/v0.1/Download/file/?uri=' +\
                         data['dataURL'][i] for i in range(len(search))]
    # data = data[data['t_exptime'].isin(cadence)]
    data = data[~data.dataURL.str.endswith('_dvt.fits')].reset_index(drop=True)
    provenance_name = data['provenance_name'].tolist()
    lc_links = data['dataURL'].tolist()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    os.makedirs('firefly/discover', exist_ok=True)
    csv_in_dir = []
    print(f'Fitting TIC {tic_id}')
    for j, fitsfile in enumerate(lc_links):
        with fits.open(fitsfile, cache=False) as TESS_fits:
            time = TESS_fits[1].data['TIME'] + 2457000
            flux = TESS_fits[1].data['PDCSAP_FLUX']
            flux_err = TESS_fits[1].data['PDCSAP_FLUX_ERR']
            
        write_dict = []
        for i in range(len(time)):
            write_dict.append({'Time': time[i], 'Flux': flux[i],
                               'Flux err': flux_err[i]})
        source = f'firefly/discover/{tic_id}'
        mast_name = data['obs_id'][j]
        os.makedirs(f'{source}', exist_ok=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
        csv_name = f'{source}/{mast_name}/{mast_name}.csv'
        write_dict = {'Time': time, 'Flux': flux, 'Flux err': flux_err}
        df = DataFrame(write_dict)
        df.to_csv(csv_name, index=False, na_rep='nan')
        csv_in_dir.append(f'{os.getcwd()}/{csv_name}')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Plot the final output
    for k, csv in enumerate(csv_in_dir):
        lc_plot(csv)

def clip(df):
    from lightkurve import LightCurve
    from astropy.stats.funcs import mad_std
    time = df['TIME']
    flux = df['PDCSAP_FLUX']
    flux_err = df['PDCSAP_FLUX_ERR']
    lc = LightCurve(time, flux, flux_err)
    lc = lc.remove_outliers(sigma_upper=6, sigma_lower=50, stdfunc=mad_std)
    df_new = DataFrame()
    df_new['TIME'] = lc.time
    df_new['PDCSAP_FLUX'] = lc.flux
    df_new['PDCSAP_FLUX_ERR'] = lc.flux_err
    print(f'Sigma clipped {len(df)-len(df_new)} cadences.')
    return df_new

def _fits(exoplanet,
          exo_folder,
          cache,
          hlsp,
          cadence,
          bitmask
):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # MAST Search
    tic_id = _tic(exoplanet).replace('TIC ', '')
    print(f'\nSearching MAST for {exoplanet} (TIC {tic_id}).')
    search = obs.query_criteria(dataproduct_type=['timeseries'],
                                project='TESS',
                                provenance_name=hlsp,
                                t_exptime=cadence,
                                target_name=tic_id
                                ).to_pandas()
    if len(search)==0:
        rmtree(exo_folder)
        print(f'Search result contains no data products for {exoplanet}.')
        sys.exit(f'Search result contains no data products for {exoplanet}.')
    data = search[['obs_id', 'target_name', 'dataURL', 't_exptime',
                   'provenance_name', 'project', 'sequence_number']]
    data['dataURL'] = ['https://mast.stsci.edu/api/v0.1/Download/file/?uri=' +\
                         data['dataURL'][i] for i in range(len(search))]
    sector_list = natsorted(data['sequence_number'].values.tolist())
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Dataframe Checks
    data = data[data['t_exptime']==cadence]
    data = data[~data.dataURL.str.endswith('_dvt.fits')].reset_index(drop=True)
    provenance_name = data['provenance_name'].tolist()
    lc_links = data['dataURL'].tolist()
    print(f'\nQuery from MAST returned {len(lc_links)} '
          f'data products for {exoplanet} (TIC {tic_id}).\n')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Print search result
    show = data[['obs_id', 'target_name', 't_exptime',
                   'provenance_name', 'project']]
    show = show.rename(columns={"obs_id": "Product",
                                "target_name": "TIC ID",
                                "t_exptime": "Cadence",
                                "provenance_name": "HLSP",
                                "project": "Mission"})
    show['Product'] = \
            Categorical(show['Product'],
            ordered=True,
            categories=natsorted(show['Product'].unique()))
    show = show.sort_values('Product')
    print(tabulate(show, tablefmt='psql', showindex=False, headers='keys'))
    print()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract Time series
    csv_in_dir = []
    for j, fitsfile in enumerate(lc_links):
        TESS_fits = None
        while TESS_fits == None:
            try:
                TESS_fits = Table.read(fitsfile, cache=cache)
            except:
                pass
            if TESS_fits != None:
                pass
        source = f'{exo_folder}/mastDownload'
        mast_name = data['obs_id'][j]
        os.makedirs(f'{source}/{mast_name}', exist_ok=True)
        TESS_fits.write(f'{source}/{mast_name}/{mast_name}.fits',
                          overwrite=True)
        TESS_fits = TESS_fits.to_pandas()
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Quality Masks
        AttitudeTweak = 1
        SafeMode = 2
        CoarsePoint = 4
        EarthPoint = 8
        Argabrightening = 16
        Desat = 32
        ApertureCosmic = 64
        ManualExclude = 128
        Discontinuity = 256
        ImpulsiveOutlier = 512
        CollateralCosmic = 1024
        #: The first stray light flag is set manually by MIT based on visual inspection.
        Straylight = 2048
        #: The second stray light flag is set automatically by Ames/SPOC based on background level thresholds.
        Straylight2 = 4096
        PlanetSearchExclude = 8192
        BadCalibrationExclude = 16384
        IsuffTargetsErrorCorr = 32768
        DEFAULT_BITMASK = [
            AttitudeTweak, SafeMode, CoarsePoint, EarthPoint, Desat, ManualExclude,
            ImpulsiveOutlier, Argabrightening, BadCalibrationExclude#, Straylight2
        ]
        HARD_BITMASK = [
            AttitudeTweak, SafeMode, CoarsePoint, EarthPoint, Desat, ManualExclude,
            ImpulsiveOutlier, Argabrightening, BadCalibrationExclude, Straylight2
        ]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
        TESS_fits['TIME'] = TESS_fits['TIME'] + 2457000
        csv_name = f'{source}/{mast_name}/{mast_name}.csv'
        _df = TESS_fits[['TIME', 'PDCSAP_FLUX', 'PDCSAP_FLUX_ERR', 'QUALITY']]
        if bitmask=='default':
            df = _df[~_df['QUALITY'].isin(DEFAULT_BITMASK)].drop('QUALITY', axis=1)
            print(f'Removed {len(_df)-len(df)} bad cadences.')
        elif bitmask=='hard':
            df = _df[~_df['QUALITY'].isin(HARD_BITMASK)].drop('QUALITY', axis=1)
            print(f'Removed {len(_df)-len(df)} bad cadences.')
            # df = clip(df)
        else:
            df = _df.drop('QUALITY', axis=1)
        df = df.rename(columns = {'TIME':'Time',
                                  'PDCSAP_FLUX':'Flux',
                                  'PDCSAP_FLUX_ERR':'Flux err'})
        df.to_csv(csv_name, index=False, na_rep='nan')
        csv_in_dir.append(f'{os.getcwd()}/{csv_name}')
    csv_in_dir = natsorted(csv_in_dir)
    print('\nSplitting up the lightcurves into seperate epochs:\n')
    # csv_in_dir = []
    # for r, d, f in os.walk(source):
    #     for item in f:
    #         if '.csv' in item:
    #             csv_in_dir.append(os.path.join(r, item))
    return csv_in_dir, sector_list
    

def _retrieval(
        # Firefly Interface
        exoplanet,
        archive='eu',
        curve_sample=1,
        clean=False,
        cache=False,
        auto=True,
        # MAST Search
        hlsp=['SPOC'],
        cadence=120,
        bitmask='default',
        # TransitFit Variables
        walks=25,
        slices=5,
        cutoff=0.25,
        window=2.5,
        nlive=300,
        fit_ttv=False,
        detrending_list=[['nth order', 2]],
        dynesty_sample='auto',
        fitting_mode='auto',
        limb_darkening_model='quadratic',
        ld_fit_method='independent',
        max_batch_parameters=25,
        batch_overlap=2,
        dlogz=None,
        maxiter=None,
        maxcall=None,
        dynesty_bounding='multi',
        normalise=True,
        detrend=True,
        detrending_limits=[[-10,10]],
        # Plotting
        plot=True,
        marker_color='dimgray',
        line_color='black',
        bin_data=True,
        binned_color='red',
):
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    exo_folder = f'firefly/{exoplanet}'
    os.makedirs(exo_folder, exist_ok=True)
    _TESS_filter()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download Archive
    if fit_ttv==True:
        archive = 'spearnet'
    if archive=='nasa':
        host_T, host_z, host_r, host_logg, t0, P, t14, repack, ra, dec, dist = \
            priors(exoplanet, archive=archive, save=True, user=False, auto=auto,
                    fit_ttv=fit_ttv)
    else:
        host_T, host_z, host_r, host_logg, t0, P, t14, repack = \
            priors(exoplanet, archive=archive, save=True, user=False, auto=auto,
                    fit_ttv=fit_ttv)
    if auto==False:
        answer = ''
        while answer!='y':
            answer = input('Modify your priors file, type y to proceed. ')
            if answer=='q':
                sys.exit('Exiting..')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Split the Light curves
    split_curve_in_dir = []
    transits_per_sector = []
    csv_in_dir, sector_list = _fits(exoplanet, exo_folder=exo_folder, cache=cache, hlsp=hlsp, cadence=cadence, bitmask=bitmask)
    for i, csvfile in enumerate(csv_in_dir):
        split_curves = split_lightcurve_file(csvfile, t0=t0, P=P, t14=t14,
                                             cutoff=cutoff, window=window)
        split_curves = [s + '.csv' for s in split_curves]
        split_curve_in_dir.append(split_curves)
        transits_per_sector.append(len(split_curves))
    split_curve_in_dir = [i for sub in split_curve_in_dir for i in sub]
    print(f'\nA total of {len(split_curve_in_dir)} lightcurves '
          'were generated.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Sort the files into ascending order and take random sample
    curves = ceil(curve_sample * len(split_curve_in_dir))
    split_curve_in_dir = random.sample(split_curve_in_dir, k=int(curves))
    split_curve_in_dir = natsorted(split_curve_in_dir)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths
    data_path = f'{exo_folder}/data_paths.csv'
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns=cols)
    for i, split_curve in enumerate(split_curve_in_dir):
        df = df.append([{'Path': split_curve}], ignore_index=True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    print(f'\nA random sample of {len(df)} lightcurves will be fitted'
          ' across all TESS Sectors.\n')
    df.to_csv(data_path, index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = data_path
    priors_csv = f'{exo_folder}/{exoplanet} Priors.csv'
    filters = 'firefly/data/TESS_filter_path.csv'
    os.makedirs('firefly/data/ldtk', exist_ok=True)
    ldtk_cache = 'firefly/data/ldtk'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'{exo_folder}/output_parameters'
    fitted_lightcurve_folder = f'{exo_folder}/fitted_lightcurves'
    plot_folder = f'{exo_folder}/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Run the retrieval
    run_retrieval(
            data,
            priors_csv,
            filters,
            host_T=host_T,
            host_logg=host_logg,
            host_z=host_z,
            host_r=host_r,
            cadence=cadence/60,
            walks=walks,
            slices=slices,
            nlive=nlive,
            fit_ttv=fit_ttv,
            detrending_list=detrending_list,
            dynesty_sample=dynesty_sample,
            fitting_mode=fitting_mode,
            limb_darkening_model=limb_darkening_model,
            ld_fit_method=ld_fit_method,
            ldtk_cache=ldtk_cache,
            max_batch_parameters=max_batch_parameters,
            batch_overlap=batch_overlap,
            dlogz=dlogz,
            maxiter=maxiter,
            maxcall=maxcall,
            dynesty_bounding=dynesty_bounding,
            normalise=normalise,
            detrend=detrend,
            detrending_limits=detrending_limits,
            results_output_folder=results_output_folder,
            final_lightcurve_folder=fitted_lightcurve_folder,
            plot_folder=plot_folder,
            plot=plot,
            marker_color=marker_color,
            line_color=line_color,
            bin_data=bin_data,
            binned_color=binned_color,
        )
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Save Best Values
    master = read_csv(f'{exo_folder}/output_parameters/Complete_results.csv',
                      index_col='Parameter')
    master = master[['Best', 'Error']]
    P = master['Best']['P']
    try:
        Perr = float(master['Error'].filter(like = 'P', axis=0))
    except ValueError:
        Perr = float()
    t0 = master['Best']['t0'][0]
    t0err = float(master['Error']['t0'][0])
    a = master['Best']['a/AU']
    aerr = float(master['Error']['a/AU'])
    ar = master['Best']['a/r*']
    arerr = float(master['Error']['a/r*'])
    rp = master['Best']['rp/r*']
    rperr = float(master['Error']['rp/r*'])
    inc = master['Best']['inc']
    incerr = float(master['Error']['inc'])
    ecc = master['Best']['ecc']
    try:
        eccerr = float(master['Error']['ecc'])
    except ValueError:
        eccerr = float()
    w = master['Best']['w']
    try:
        werr = float(master['Error']['w'])
    except ValueError:
        werr = float()
    try:
        mad, madbin = density_scatter(exoplanet=exoplanet, transits=int(len(df)))
    except Exception as e:
        print(e)
    t_depth = rp**2
    sens = t_depth/mad
    sens_bin = t_depth/madbin
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    print(
        'Variables used:\n\n'
        f'target={exoplanet}\n'
        f'archive={archive}\n'
        f'hlsp={hlsp}\n'
        f'curve_sample={str(curve_sample)}\n'
        f'clean={clean}\n'
        f'cache={cache}\n'
        f'cutoff={str(cutoff)}\n'
        f'window={str(window)}\n'
        f'nlive={str(nlive)}\n'
        f'fit_ttv={fit_ttv}\n'
        f'detrending_list={str(detrending_list)}\n'
        f'dynesty_sample={dynesty_sample}\n'
        f'fitting_mode={fitting_mode}\n'
        f'limb_darkening_model={limb_darkening_model}\n'
        f'ld_fit_method={ld_fit_method}\n'
        f'max_batch_parameters={str(max_batch_parameters)}\n'
        f'batch_overlap={str(batch_overlap)}\n'
        f'dlogz={str(dlogz)}\n'
        f'maxiter={str(maxiter)}\n'
        f'maxcall={str(maxcall)}\n'
        f'dynesty_bounding={dynesty_bounding}\n'
        f'normalise={normalise}\n'
        f'detrend={detrend}',
        file=open(exo_folder+'/variables.txt', 'w')
    )
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    sci_prod = ' '.join(hlsp)
    archive_name = f"{exoplanet} {archive.upper()} {sci_prod} {now}"
    if clean==True:
        try:
            rmtree(f'{exo_folder}/output_parameters/quicksaves')
            rmtree(f'{exo_folder}/output_parameters/filter_0_parameters/quicksaves')
        except Exception:
            pass
    if fit_ttv==True:
        try:
            here = os.path.dirname(os.path.abspath(__file__))
            spearnet_csv = f'{here}/data/spear.csv'
            spearnet = read_csv(spearnet_csv).set_index('pl_name')
            s = spearnet.loc[[exoplanet]]
            t0ttv = s['pl_tranmid'][0]
            t0errttv = s['pl_tranmiderr1'][0]
            file = f'firefly/{exoplanet}/output_parameters/Complete_results.csv'
            chi2_red, nsig, loss, fap, period = oc_fold(t0ttv, t0errttv,
                                                        file=file, exoplanet=exoplanet,
                                                        transits_per_sector=transits_per_sector,
                                                        sector_list=sector_list)
            data = {'pl_name':exoplanet, 'red_chi2':chi2_red, 'sigma':nsig,
                    'mean_avg_err':loss,
                    'fap':fap, 'o-c_period':period,
                    'pl_orbper':P, 'pl_orbpererr1':Perr,
                    'pl_tranmid':t0, 'pl_tranmiderr1':t0err,
                    'pl_orbsmax':a, 'pl_orbsmaxerr1':aerr, 'pl_radj':rp,
                    'pl_radjerr1':rperr, 'pl_orbincl':inc,
                    'pl_orbinclerr1':incerr, 'pl_orbeccen':ecc, 'pl_orbeccenerr1':eccerr,
                    'pl_orblper':w, 'pl_orblpererr1':werr, 'Transits':int(len(df)),
                    'Date':now, 'Archive':archive.upper(), 'Unbinned Sigma':mad,
                    'Binned Sigma':madbin,
                'Transit Depth':t_depth, 'Sensitivity':sens,
                'Binned Sensitivity':sens_bin
            }
            df = DataFrame(data, index=[0])
            summary_master = 'firefly/data/spear_ttv.csv'
            if not os.path.exists(summary_master):
                df.to_csv(summary_master, index=False)
            else:
                add = read_csv(summary_master)
                add = add.append(df)
                add['pl_name'] = \
                    Categorical(add['pl_name'],
                    ordered=True,
                    categories=natsorted(add['pl_name'].unique()))
                add = add.sort_values('pl_name')
                add .to_csv(summary_master, index=False)
        except Exception as e:
            print(e)
        os.makedirs('firefly/ttv', exist_ok=True)
        os.makedirs(f'firefly/ttv/{fitting_mode}', exist_ok=True)
        archive_folder = f'firefly/ttv/{fitting_mode}/{archive_name}'
        make_archive(archive_folder, format='zip',
                     root_dir=f'{os.getcwd()}/firefly/',
                     base_dir=f'{exoplanet}')
    else:
        data = {'pl_name':exoplanet, 'pl_orbper':P, 'pl_orbpererr1':Perr,
                'pl_tranmid':t0, 'pl_tranmiderr1':t0err,
                'pl_orbsmax':a, 'pl_orbsmaxerr1':aerr, 'pl_radj':rp,
                'pl_radjerr1':rperr, 'pl_orbincl':inc,
                'pl_orbinclerr1':incerr, 'pl_orbeccen':ecc, 'pl_orbeccenerr1':eccerr,
                'pl_orblper':w, 'pl_orblpererr1':werr, 'Transits':int(len(df)),
                'Date':now, 'Archive':archive.upper(), 'Unbinned Sigma':mad,
                'Binned Sigma':madbin,
                'Transit Depth':t_depth, 'Sensitivity':sens,
                'Binned Sensitivity':sens_bin
        }
        df = DataFrame(data, index=[0])
        summary_master = 'firefly/data/spear.csv'
        if not os.path.exists(summary_master):
            df.to_csv(summary_master, index=False)
        else:
            add = read_csv(summary_master)
            add = add.append(df)
            add['pl_name'] = \
                Categorical(add['pl_name'],
                ordered=True,
                categories=natsorted(add['pl_name'].unique()))
            add = add.sort_values('pl_name')
            add .to_csv(summary_master, index=False)
        os.makedirs(f'firefly/{fitting_mode}', exist_ok=True)
        archive_folder = f'firefly/{fitting_mode}/{archive_name}'
        make_archive(archive_folder, format='zip',
                     root_dir=f'{os.getcwd()}/firefly/',
                     base_dir=f'{exoplanet}')
    rmtree(f'{exo_folder}')
    return archive_name


