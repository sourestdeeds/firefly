# retrieval

A target data retriever for confirmed/candidate TESS exoplanets.
Generates the priors and host star variables for a chosen target.
Downloads exoplanet archives every 10 days and stores in /data.
Target lightcurve files are downloaded from MAST, fits file is
stored in Planet/exoplanet/exoplanet.fits.

## Dependancies
```
lightkurve
transitfit
```
An example use with TransitFit is the following:
```
exoplanet = 'WASP-43 b'
retrieval(exoplanet)
```
# auto_retrieval

Automated version of retrieval. Sends an email to transitfit.server@gmail.com
upon an error or full completion of a target. Iteratively takes targets and
employs TransitFit across each TESS sector for every exoplanet in the list given.
Runs TransitFit for all available split curves.

An example use with TransitFit is the following:
```
file = ('WASP-43 b', 'WASP-18 b', 'WASP-91 b', 'WASP-12 b',
            'WASP-126 b', 'LHS 3844 b', 'GJ 1252 b', 'TOI-270 b')
auto_retrieval(file)
```

    Parameters
    ----------
    file : str
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

    Returns
    -------
    A whole lot of data to science!
