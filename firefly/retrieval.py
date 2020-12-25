#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A target data retriever for confirmed/candidate TESS exoplanets.

@author: Steven Charles-Mindoza
"""

from ._utils import _fits, _TESS_filter, _MAST_query, _email
from ._archive import _eu, _nasa, _check_nan
from ._search import _fuzzy_search
from .query import tess

from transitfit import split_lightcurve_file, run_retrieval
try:
    from lightkurve import search_lightcurve
except:
    try:
        from lightkurve import search_lightcurvefile as search_lightcurve
    except:
        raise
from traceback import format_exc
from tabulate import tabulate
from datetime import datetime
from pandas import DataFrame
from shutil import rmtree, make_archive
import sys
import os
import random


def _nan(exoplanet, archive, printing=False):
    nan = _check_nan(exoplanet, archive=archive)
    if nan == True:
        _check_nan(exoplanet, archive=archive, printing=True)
        verify = ''
        while (verify!="y" and verify!="n"):
            verify = input(f'\nWARNING: {exoplanet} has missing '
                                'prior entries. Proceed ([y]/n)?\n')
        if verify == "n":
            sys.exit()
        elif verify == "y":
            pass


def _retrieval_input_target(exoplanet, archive):
    # Check inputs are sensible
    if not (archive == 'eu' or archive == 'nasa'):
        sys.exit('Archive data options are: \'eu\' or \'nasa\'')
    archive_data = 'EU archive'
    highest, ratios = _fuzzy_search(exoplanet, archive=archive)
    exoplanet = highest[0]
    verify = ''
    while (verify!="y" and verify!="n" and verify!='q'):
        print(f'\nTarget search chose {exoplanet}.')
        verify = input('Proceed ([y]/n)?\n')
        _nan(exoplanet, archive)
    if (verify=='q'):
        sys.exit('You chose to quit.')
    elif (verify=='y'):
        print(f'\nChecking data products from MAST for {exoplanet}.')
        return highest[0]
    elif (verify=='n'):
        while (verify!="y" and verify!='q' and exoplanet!='q'):
            tess()
            exoplanet = input('Please refine your search or type q to quit: ')
            highest, ratios = _fuzzy_search(exoplanet, archive=archive)
            exoplanet = highest[0]
            print(f'\nTarget search chose {highest[0]} from the '
                  f'{archive_data}:\n')
            verify = input('Proceed ([y]/n)? or type q to quit.\n')
            _nan(exoplanet, archive)
        if (verify=='q' or exoplanet=='q'):
            sys.exit('You chose to quit.')
        elif (verify=="y"):
            print(f'\nChecking data products from MAST for {exoplanet}.')
            return highest[0]


def retrieval_input_sector(sector_list):
    sector = '0'
    check = all(i in sector_list for i in sector) # False
    while (check!=True and sector!='q' and sector!='a'):
        print('\nAvailable TESS Sectors are:', *sector_list)
        sector = input('Enter which TESS Sectors you would like to download' 
                       ' or type a to download them all: ')
        check = all(item in sector_list for item in sector.split())
    if sector == 'q':
            sys.exit('You chose to quit.')
    elif sector == 'a':
        return sector_list
    elif check == True:
        sector_list = sector
        return sector_list


def retrieval_input_curve_sample(split_curve_in_dir):
    sample = '0'
    sample_range = range(1, len(split_curve_in_dir)+1)
    sample_range = [str(sample) for sample in sample_range]
    while (sample not in sample_range and sample != 'q'):
        print(f'\nUp to {len(split_curve_in_dir)} epochs are available.')
        sample = input('Enter a random sample of lightcurves you wish to fit: ')
    if sample == 'q':
        sys.exit('You chose to quit.')
    elif sample in sample_range:
        return sample


def retrieval(exoplanet, archive='eu', email=False, 
              to=['transitfit.server@gmail.com'], clean=False, nlive=300, 
              fit_ttv=False, detrending_list=[['nth order', 2]],
              dynesty_sample='auto', fitting_mode='folded',
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
        >>> retrieval(exoplanet)
        
    Input is capable of handling:
    
        >>> 'wasp43b', 'WASp43b' etc
        
    Forces corrections based on classifier: 
    
        >>> 'WASP', 'LTT', 'GJ' etc

    Parameters
    ----------
    exoplanet : 'WASP-43 b', string, 
        The target exoplanet.

    archive : string, optional
        Defines which exoplanet archive the priors are generated from.
        Allows for inputs 'nasa' or 'eu'. The default is 'nasa'.

        EU : Data downloaded and stored in 'data/eu.csv'.
        http://exoplanet.eu/catalog/#

        NASA : Data downloaded and stored in 'data/nasa.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html
        
    email : bool, optional
        If True will send status emails. The default is False.
    to : str, optional
        The email address to send status updates to.
        The default is:
            
        >>> to=['transitfit.server@gmail.com']
    clean : bool, optional
        If True will delete all downloaded files and zip outputs only.
        The default is False.
    nlive : int, optional
        The number of live points to use in the nested sampling retrieval.
        Default is 300.
    detrending_list : array_like, shape (n_detrending_models, 2)
        A list of different detrending models. Each entry should consist
        of a method and a second parameter dependent on the method.
        Accepted methods are:
            
        - ['nth order', order]
        - ['custom', function, [global fit indices, filter fit indices, 
                                epoch fit indices]]
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
    if not (archive == 'eu' or archive == 'nasa'):
        sys.exit('Archive data options for dtype are: \'eu\' or \'nasa\'')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Target Input
    exoplanet = _retrieval_input_target(exoplanet, archive)
    exo_folder = f'firefly/{exoplanet}'
    try:
        rmtree(exo_folder)
    except:
        pass
    os.makedirs(exo_folder, exist_ok=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Filter Setup
    _TESS_filter()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Search MAST lightcurves
    sector_list = _MAST_query(exoplanet, exo_folder)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Sector Input
    sector_list = retrieval_input_sector(sector_list)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Download Archive
    if archive == 'eu':
        host_T, host_z, host_r, host_logg, t0, P, t14, nan = _eu(exoplanet)
    elif archive == 'nasa':
        host_T, host_z, host_r, host_logg, t0, P, t14, nan = _nasa(exoplanet)
    cols = [['t0', t0], ['P', P], ['t14', t14]]
    df = DataFrame(cols, columns=['Parameter', 'Value'])
    print('\nSplitting the lightcurve into seperate epochs'
          ' using the following parameters.\n')
    print(tabulate(df, tablefmt='psql', showindex=False, headers='keys'))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Iterate over all sectors
    for i, sector in enumerate(sector_list):
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # MAST Download
        lc = search_lightcurve(exoplanet, mission='TESS', #radius=750,
                                   sector=sector)
        print(f'\nDownloading MAST Lightcurve for {exoplanet} -' +
              f' TESS Sector {str(sector)}.')
        lc.download_all(download_dir=exo_folder)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Split the Light curves
    split_curve_in_dir = []
    csv_in_dir = _fits(exoplanet, exo_folder, clean)
    for i, csvfile in enumerate(csv_in_dir):
        split_curves = split_lightcurve_file(csvfile, t0=t0, P=P) #, t14=t14)
        split_curves = [s + '.csv' for s in split_curves]
        split_curve_in_dir.append(split_curves)
        if clean == True:
            os.remove(csvfile)
    split_curve_in_dir = [i for sub in split_curve_in_dir for i in sub]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Sort the files into ascending order and take a random sample
    curve_sample = retrieval_input_curve_sample(split_curve_in_dir)
    split_curve_in_dir = random.sample(split_curve_in_dir, k=int(curve_sample))
    split_curve_in_dir.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Set the data paths
    data_path = f'{exo_folder}/data_paths.csv'
    cols = ['Path', 'Telescope', 'Filter', 'Epochs', 'Detrending']
    df = DataFrame(columns=cols)
    for i, split_curve in enumerate(split_curve_in_dir):
        df = df.append([{'Path': split_curve}], ignore_index=True)
        df['Telescope'], df['Filter'], df['Detrending'] = 0, 0, 0
        df['Epochs'] = range(0, len(df))
    print(f'\nIn total, {len(df)} lightcurves will be fitted across all'
          ' TESS Sectors.\n')
    df.to_csv(data_path, index=False, header=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Paths to data, priors, and filter info:
    data = data_path
    priors = f'{exo_folder}/{exoplanet} Priors.csv'
    filters = 'firefly/data/TESS_filter_path.csv'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Output folders
    results_output_folder = f'{exo_folder}/output_parameters'
    fitted_lightcurve_folder = f'{exo_folder}/fitted_lightcurves'
    plot_folder = f'{exo_folder}/plots'
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Run the retrieval
    try:
        run_retrieval(data, priors, filters, 
                      detrending_list=detrending_list,
                      host_T=host_T, host_logg=host_logg, 
                      host_z=host_z, host_r=host_r, 
                      dynesty_sample=dynesty_sample,
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
    except KeyboardInterrupt:
        now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
        keyboard = f'firefly/KeyboardInterrupt/{exoplanet} ' +\
                   f'{now} KeyboardInterrupt'
        make_archive(keyboard, format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{exoplanet}')
        rmtree(exo_folder)
        raise
    except BaseException:
        now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
        exception = f'firefly/Exception/{exoplanet} ' +\
                    f'{now} Exception'
        make_archive(exception, format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{exoplanet}')
        rmtree(exo_folder)
        trace_back = format_exc()
        if email == True:
            _email(f'Exception: {exoplanet}', trace_back, to=to)
        raise
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Cleanup
    if clean == True:
        rmtree(f'{exo_folder}/mastDownload')
        os.remove(data)
        os.remove(priors)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Archive and sort
    try:
        now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
        make_archive(f'{exo_folder} {now}', format='gztar',
                     root_dir=f'{os.getcwd()}/firefly/',
                     base_dir=f'{exoplanet}')
        rmtree(exo_folder)
        success = f'{os.getcwd()}/{exo_folder} {now}.gz.tar'
        print(f'\nData location: {success}\n'
                           'A new target has been fully retrieved across ' +
                           'all available TESS Sectors.')
        if email == True:
            _email(f'Success: {exoplanet}',
                           f'Data location: {success} \n\n'
                           'A new target has been fully retrieved across ' +
                           'all available TESS Sectors.', to=to)
    except:
        pass
