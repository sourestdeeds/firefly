#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A target data retriever for confirmed/candidate TESS exoplanets.

@author: Steven Charles-Mindoza
"""

from .query import query, tess_targets
from ._utils import _fits, _TESS_filter
from ._archive import _eu, _nasa
from ._search import _fuzzy_search

from transitfit import split_lightcurve_file, run_retrieval
try:
    from lightkurve import search_lightcurve
except:
    try:
        from lightkurve import search_lightcurvefile as search_lightcurve
    except:
        raise
from tabulate import tabulate
from datetime import datetime
from pandas import DataFrame
from shutil import rmtree, move, make_archive
import sys
import os


def _retrieval_input_target(target, archive):
    # Check inputs are sensible
    if not (archive == 'eu' or archive == 'nasa'):
        sys.exit('Archive data options are: \'eu\' or \'nasa\'')
    archive_data = 'EU archive'
    highest, ratios = _fuzzy_search(target, archive=archive)
    target = highest[0]
    verify = ''
    while (verify!="y" and verify!="n" and verify!='q' and verify!='tess'):
        print(f'\nTarget search chose {target} from the {archive_data}:\n')
        query(target, archive=archive)
        print('\nFor a list of TESS targets, type tess.')
        verify = input('Proceed ([y]/n)?\n')
    if (verify=='q'):
        sys.exit('You chose to quit.')
    elif (verify=='y'):
        print(f'\nChecking data products from MAST for {target}.')
        pass
        return highest[0]
    elif (verify =='tess' or verify=='n'):
        while (verify!="y" and verify!='q'):
            tess_targets()
            target = input('Please refine your search: ')
            highest, ratios = _fuzzy_search(target, archive=archive)
            target = highest[0]
            print(f'\nTarget search chose {highest[0]} from the '
                  f'{archive_data}:\n')
            query(target, archive=archive)
            verify = input('Proceed ([y]/n)? or type q to quit.\n')   
        if (verify=='q'):
            sys.exit('You chose to quit.')
        elif (verify=="y"):
            print(f'\nChecking data products from MAST for {target}.')
            pass
            return highest[0]


def retrieval(target, archive='eu', nlive=300, fit_ttv=False,
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
    Target lightcurve files are downloaded from MAST, then split into 
    separate epochs. Upon user entry of the amount of epochs to fit,
    TransitFit will fit the curves and return the results. The results
    are then zipped up and time stamped.

    An example use with TransitFit is the following:

        >>> from firefly import retrieval
        >>> exoplanet = 'WASP-43 b'
        >>> retrieval(target)
        
    Input is capable of handling:
    
        >>> 'wasp43b', 'WASp43b' etc
        
    Forces corrections based on classifier: 
    
        >>> 'WASP', 'LTT', 'GJ' etc

    Parameters
    ----------
    target : 'WASP-43 b', string, 
        The target exoplanet.

    archive : string, optional
        Defines which exoplanet archive the priors are generated from.
        Allows for inputs 'nasa' or 'eu'. The default is 'nasa'.

        EU : Data downloaded and stored in 'data/eu.csv'.
        http://exoplanet.eu/catalog/#

        NASA : Data downloaded and stored in 'data/nasa.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html

    nlive : int, optional
        The number of live points to use in the nested sampling retrieval.
        Default is 300.
    detrending_list : array_like, shape (n_detrending_models, 2)
        A list of different detrending models. Each entry should consist
        of a method and a second parameter dependent on the method.
        Accepted methods are:
            
        - ['nth order', order]
        - ['custom', function, [global fit indices, filter fit indices, epoch fit indices]]
        - ['off', ]
        Function here is a custom detrending function. TransitFit assumes
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
            
        >>> foo(times, a, b, c):
                # do something
           
        and a should be fitted globally, then the entry in the method_list
        would be 
        
        - ['custom', foo, [1], [], []].
    dynesty_sample : str, optional
        Method used to sample uniformly within the likelihood constraint,
        conditioned on the provided bounds. Unique methods available are:
            
        - uniform sampling within the bounds('unif') 
        - random walks with fixed proposals ('rwalk') 
        - random walks with variable ("staggering") proposals ('rstagger') 
        - multivariate slice sampling along preferred orientations ('slice') 
        - "random" slice sampling along all orientations ('rslice') 
        - "Hamiltonian" slices along random trajectories ('hslice') 
        and any callable function which follows
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
        
        >>> ln(z + z_est) - ln(z) < dlogz, 
        
        where z is the current evidence
        from all saved samples and z_est is the estimated contribution from
        the remaining volume. The default is 
        
        >>> 1e-3 * (nlive - 1) + 0.01.
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
            
        >>> 1/f_min <= c_n <= 1/f_max
        
        as the default range, where f_min and f_max are the minimum and maximum
        flux values for a given light curve. Default is True.
    detrend : bool, optional
        If True, will initialise detrending fitting. Default is True.

    Returns
    -------
    Zipped files are found in :
    
    >>> firefly/WASP-43 b timestamp.gz.tar
    '''
    target = _retrieval_input_target(target, archive)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    _TESS_filter()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download MAST lightcurves
    lc = search_lightcurve(target, mission='TESS')
    if len(lc) == 0:
        sys.exit(f'Search result contains no data products for {target}.')
    sector_list = lc .table .to_pandas()['sequence_number'] \
                     .drop_duplicates() .tolist()
    sector_list = [str(sector) for sector in sector_list]
    print(f'\nQuery from MAST returned {len(sector_list)} '
          f'data products for {target}.\n')
    lc = lc .table .to_pandas()[['observation', 
                                 'productFilename', 'size', 't_exptime']] \
            .rename(columns={'observation':'Observation'}) \
            .rename(columns={'size':'Size'}) \
            .rename(columns={'productFilename':'Product'}) 
    lc = lc[lc.t_exptime != 20].drop(['t_exptime'], axis=1)
    print(tabulate(lc, tablefmt='psql', showindex=False, headers='keys'))
    sector = '0'
    while (sector not in sector_list and sector != 'q'):
        print('\nAvailable TESS Sectors are:', *sector_list)
        sector = input('Enter which TESS Sector you would like to download: ')
        while (sector not in sector_list and sector != 'q'):
            print('\nPlease choose an integer from the list: \n', *sector_list,
                  '\nor type q to quit.')
            sector = input('Enter which TESS Sector you' +
                               ' would like to download: ')
        if sector == 'q':
            sys.exit('You chose to quit.')
        elif sector in sector_list:
            pass
    sector = int(sector)
    exo_folder = f'firefly/{target}/TESS Sector {str(sector)}'
    os.makedirs(exo_folder, exist_ok=True)
    lc = search_lightcurve(target, mission='TESS',
                               sector=sector)
    print(f'\nDownloading MAST Lightcurve for TESS Sector {str(sector)}.')
    lc.download_all(download_dir=exo_folder)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Extract all light curves to a single csv file
    fitsfile = _fits(target, exo_folder)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download Archive
    if archive == 'eu':
        host_T, host_z, host_r, host_logg, t0, P, t14, nan = _eu(target)
    elif archive == 'nasa':
        host_T, host_z, host_r, host_logg, t0, P, t14, nan = _nasa(target)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Split the Light curves
    cols = [['t0', t0], ['P', P], ['t14', t14]]
    df = DataFrame(cols, columns=['Parameter', 'Value'])
    print('\nSplitting the lightcurve into seperate epochs'
          ' using the following parameters.\n')
    print(tabulate(df, tablefmt='psql', showindex=False, headers='keys'))
    csvfile = f'{exo_folder}/{target}.csv'
    split_curves = split_lightcurve_file(csvfile, t0=t0, P=P, t14=t14)
    print(f'\nA total of {str(len(split_curves))} lightcurves were created.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the Data Paths
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns=cols)
    curves = '0'
    curve_list = [str(curve) for curve in range(1,len(split_curves)+1)]
    while (curves not in curve_list and curves != 'q'):
        curves = input('Enter how many lightcurves you wish to fit: ')
        while (curves not in curve_list and curves != 'q'):
            print('\nPlease enter an integer between 1 and' +\
                  f' {len(split_curves)} or type q to quit.')
            curves = input('Enter how many lightcurves you wish to fit: ')
        if curves == 'q':
            sys.exit('You chose to quit.')
        elif (curves in curve_list):
            pass
    print()
    for i in range(int(curves)):
        df = df.append([{'Path': f'{os.getcwd()}/{exo_folder}' +
                         f'/split_curve_{str(i)}.csv'}],
                       ignore_index=True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    df.to_csv(r'firefly/data/data_paths.csv', index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = 'firefly/data/data_paths.csv'
    priors = f'firefly/{target}/{target} Priors.csv'
    filters = 'firefly/data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'{exo_folder}/output_parameters'
    fitted_lightcurve_folder = f'{exo_folder}/fitted_lightcurves'
    plot_folder = f'{exo_folder}/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Cleanup
    try:
        rmtree(f'{exo_folder}/mastDownload')
        move(fitsfile, f'{exo_folder}.fits')
        os.remove(f'{exo_folder}.fits')
        os.remove(csvfile)
        os.remove(priors)
        for i in range(len(split_curves)):
            os.remove(f'{exo_folder}/split_curve_{str(i)}.csv')
    except BaseException:
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
    make_archive(f'{exo_folder} {now}', format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{target}')
    rmtree(f'{exo_folder}')
    success = f'{os.getcwd()}/{exo_folder} {now}.gz.tar'
    print(f'\nData location: {success}\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.')
