#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated version of retrieval.

@author: Steven Charles-Mindoza
"""

from ._utils import _email, _retrieval
from ._archive import _search, _load_csv

from datetime import datetime
from shutil import rmtree, make_archive
from traceback import format_exc
import sys
import os



def _auto_input_check(targets, curve_sample):
    if not (0 < curve_sample <= 1):
        sys.exit('The curve sample must be in the range 0 < curve_sample <= 1.')
    _load_csv()
    targets = ''.join(targets)
    highest, ratios = _search(targets)
    exoplanet = highest[0]
    print(f'Target search chose {exoplanet}.')
    return exoplanet


def firefly(
        # Firefly Interface
        targets,
        archive='eu',
        hlsp=['SPOC'],
        curve_sample=1,
        email=False,
        to=['transitfit.server@gmail.com'],
        clean=True,
        cache=False,
        auto=True,
        # TransitFit Variables
        cutoff=0.25,
        window=2.5,
        nlive=1000,
        fit_ttv=False,
        detrending_list=[['nth order', 2]],
        dynesty_sample='rslice',
        fitting_mode='folded',
        limb_darkening_model='quadratic',
        ld_fit_method='coupled',
        max_batch_parameters=25,
        batch_overlap=2,
        dlogz=None,
        maxiter=None,
        maxcall=None,
        dynesty_bounding='multi',
        normalise=True,
        detrend=True
):
    '''
    Automated version of retrieval. For a single target the procedure is:
         
         >>> from firefly import firefly
             target = ('WASP-43 b',)
             firefly(target)
         
    For a list of targets:
         
         >>> from firefly import firefly
             targets = ('WASP-43 b', 'WASP-18 b')
             firefly(targets)
     
    - Targets passed are corrected for basic user input; 'wasp43b' is
      interpreted as 'WASP-43 b'. List must be of the form given in the
      example below.
    - Initial checks for targets from the exoplanet archive are then taken
      to ascertain
      whether the prior data extracted has entries in all columns. If there
      are missing entries for a given target in the list, the user will be
      asked whether to proceed.
    - Iteratively takes the targets given and employs TransitFit across each
      TESS sector for every exoplanet in the list given.
    - All available split curves are fitted with TransitFit, then the results
      are then zipped up and time stamped. Optionally sends an email upon an
      error or full completion of a target.
    - Email is also disabled by default. If enabled, status updates on
      completion and exceptions with the full traceback are sent.
      
    Background tasks for feeding data to TransitFit include:
    
    - Set the filter path to the TESS Filter.
    - Download EU/NASA exoplanet archive data every 10 days
      (checks the file age).
    - Download MAST lightcurve data for target TESS Sectors.
    - Split the lightcurves into seperate transits or epochs.
    - Create the data paths to each seperate epoch.
    - Run TransitFit.
    - Delete all downloaded data and inputs.
    - Zip and timestamp the output.
        
    Input is capable of handling :
    
        >>> ('wasp43b', 'WASp18b', 'wasP91--b') etc
        
    Forces corrections based on classifier:
    
        >>> 'WASP', 'LTT', 'GJ' etc
    
    Parameters
    ----------
    targets : str, list
        A list of exoplanet targets.
        Input is a list tuple of strings:
            
        >>> ('WASP-43 b', 'WASP-18 b', 'WASP-91 b')
    curve_sample : int {0 < curve_sample <= 1}, optional
        The fraction of curves generated to fit against. For example, setting
        
        >>> curve_sample = 0.5
        
        will fit half the curves extracted. The formula for this works as:
            
        >>> total_curves =
            curves_extracted * curve_sample
        
        Always returns an int. For example:
            
        >>> curve_sample = 0.001
        
        will fit using only 1 lightcurve from each sector.
        The default is 1 to fit all lightcurves across all sectors.
    archive : str, {'eu', 'nasa', 'org', 'all'}
        The archive to generate priors from. All takes the IQR of all
        archives (including OEC) and then the mean.
        The default is 'eu'.
    hlsp : str list, ['SPOC', 'TESS-SPOC', 'TASOC']
        SPOC is the primary TESS mission, and the rest are HLSP.
        The default is ['SPOC'].
    email : bool, optional
        If True will send status emails. The default is False.
    to : str, optional
        The email address to send status updates to.
        
        >>> to=['transitfit.server@gmail.com']
    clean : bool, optional
        If True will delete all downloaded files and zip outputs only.
        The default is False.
    cache : bool, optional
        If True will cache lightcurve fits files in root/.astropy.
        The default is False.
    auto : bool, optional
        If False will allow the user to modify the generated priors file
        before proceeding.
        The default is True.
    cutoff : float, optional
        If there are no data within
        
        >>> t14 * cutoff of t0,
        
        a period will be
        discarded. Default is 0.25
    window : float, optional
        Data outside of the range
        
        >>> [t0 ± (0.5 * t14) * window]
        
        will be discarded.
    nlive : int, optional
        The number of live points to use in the nested sampling retrieval.
        Default is 1000.
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
        the problem (from ndim).
        - When ndim < 10, this defaults to 'unif'.
        - When 10 <= ndim <= 20, this defaults to 'rwalk'.
        - When ndim > 20, this defaults to 'hslice' if a gradient is provided
          and 'slice' otherwise.
        - 'rstagger' and 'rslice' are provided as alternatives for
          'rwalk' and 'slice', respectively. Default is 'auto'.
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
        - `'2_stage'` : Fits in 2 stages, first detrending the light curves and
          then fitting the detrended curves simultaneously, using the
          `'batched'` approach if required.
        - `'batched'` : Useful for large numbers of light curves with
          relatively few shared filters, so `'folded'` loses large amounts of
          multi-epoch information. This mode splits the filters into sets of
          overlapping batches, runs each batch and uses the weighted means of
          each batch to produce a final result.
        Default is `'auto'`.
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
    A whole lot of data to science!
    Zipped files are found in:
    
    >>> firefly/WASP-43 b timestamp.gz.tar

    '''
    exoplanet = _auto_input_check(targets, curve_sample=curve_sample)
    try:
        archive_name, repack, results = \
        _retrieval(
            # Firefly Interface
            exoplanet,
            archive=archive,
            hlsp=hlsp,
            curve_sample=curve_sample,
            clean=clean,
            cache=cache,
            auto=auto,
            # TransitFit Variables
            cutoff=cutoff,
            window=window,
            nlive=nlive,
            fit_ttv=fit_ttv,
            detrending_list=detrending_list,
            dynesty_sample=dynesty_sample,
            fitting_mode=fitting_mode,
            limb_darkening_model=limb_darkening_model,
            ld_fit_method=ld_fit_method,
            max_batch_parameters=max_batch_parameters,
            batch_overlap=batch_overlap,
            dlogz=dlogz,
            maxiter=maxiter,
            maxcall=maxcall,
            dynesty_bounding=dynesty_bounding,
            normalise=normalise,
            detrend=detrend
        )
        success = f'{os.getcwd()}/firefly/{archive_name}.zip'
        print(f'\nData location: {success}\n'
                'A new target has been fully retrieved across ' +
                'all available TESS Sectors.')
        if email == True:
            exo_folder = f'firefly/{exoplanet}'
            _email(
            f'Success: {exoplanet}',
            f'Data location: {success} <br><br>'
            f'Priors from the Archive for {exoplanet}:<br>' +
            repack.to_html(float_format=lambda x: '%10.5f' % x) + '<br>'
            'TransitFit Complete Results:<br>' +
            results.to_html(float_format=lambda x: '%10.5f' % x) + '<br>'
            'Variables used:<br><br>'
            f'target={exoplanet}<br>'
            f'archive={str(archive)}<br>'
            f'hlsp={hlsp}<br>'
            f'curve_sample={str(curve_sample)}<br>'
            f'clean={clean}<br>'
            f'cache={cache}<br>'
            f'cutoff={str(cutoff)}<br>'
            f'window={str(window)}<br>'
            f'nlive={str(nlive)}<br>'
            f'fit_ttv={fit_ttv}<br>'
            f'detrending_list={str(detrending_list)}<br>'
            f'dynesty_sample={dynesty_sample}<br>'
            f'fitting_mode={fitting_mode}<br>'
            f'limb_darkening_model={limb_darkening_model}<br>'
            f'ld_fit_method={ld_fit_method}<br>'
            f'max_batch_parameters={str(max_batch_parameters)}<br>'
            f'batch_overlap={str(batch_overlap)}<br>'
            f'dlogz={str(dlogz)}<br>'
            f'maxiter={str(maxiter)}<br>'
            f'maxcall={str(maxcall)}<br>'
            f'dynesty_bounding={dynesty_bounding}<br>'
            f'normalise={normalise}<br>'
            f'detrend={detrend}',
            to=to)
    except KeyboardInterrupt:
        exo_folder = f'firefly/{exoplanet}'
        now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
        keyboard = f'firefly/KeyboardInterrupt/{exoplanet} ' +\
                   f'{now} KeyboardInterrupt'
        make_archive(keyboard, format='gztar',
                 root_dir=f'{os.getcwd()}/firefly/',
                 base_dir=f'{exoplanet}')
        rmtree(exo_folder)
        raise
        sys.exit()
    except BaseException:
        try:
            exo_folder = f'firefly/{exoplanet}'
            now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
            exception = f'firefly/Exception/{exoplanet} ' +\
                        f'{now} Exception'
            trace_back = format_exc()
            print(trace_back, file=open(exo_folder+'/traceback.txt', 'w'))
            print(
                'Variables used:\n\n'
                f'target={exoplanet}\n'
                f'archive={str(archive)}\n'
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
            make_archive(exception, format='gztar',
                     root_dir=f'{os.getcwd()}/firefly/',
                     base_dir=f'{exoplanet}')
            rmtree(exo_folder)
            if email == True:
                _email(f'Exception: {exoplanet}', trace_back, to=to)
        except:
            pass
