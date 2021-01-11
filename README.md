# **firefly**
A target selector for use with TransitFit to fit TESS lightcurves.

#### Installation (add later)
```bash
python setup.py install
```

#### Dependancies
```
transitfit
tabulate
python-Levenshtein
fuzzywuzzy
numpy
pandas
```

- Targets passed are corrected for basic user input. 'wasp43b' is
interpreted as 'WASP-43 b'. List must be of the form given in the example below.
- Initial checks for targets from the exoplanet archive are then taken to ascertain 
whether the prior data extracted has entries in all columns. If there are missing
entries for a given target in the list, the user will be asked whether to proceed.
- Iteratively takes the targets given and employs TransitFit across each TESS sector 
for every exoplanet in the list given.
- All available split curves are fitted with TransitFit, then the results
are then zipped up and time stamped. Optionally sends an email upon an error or 
full completion of a target.
- The printing of all but essential output such as exceptions are disabled. 
If you wish to run the fitting procedure for a list in piecemeal across a 
single core, set the processess variable to 1 and the list will be worked 
through one by one. It is advised that you only set the variable printing 
to True when fitting in this manner. Multiple cores give a chaotic output, 
and impacts performance.
- Email is also disabled by default. If enabled, status updates on completion
and exceptions with the full traceback are sent. (Unfinished)

Background tasks for feeding data to TransitFit include:
- Set the filter path to the TESS Filter.
- Download EU/NASA exoplanet archive data every 10 days (checks the file age).
- Download MAST lightcurve data for target TESS Sectors.
- Split the lightcurves into seperate transits or epochs.
- Create the data paths to each seperate epoch.
- Run TransitFit.
- Delete all downloaded data and inputs.
- Zip up and timestamp the output.

An example use with TransitFit is the following:
```python
from firefly import firefly, tess_viable

# For a single target
target = ('WASP-43 b',)
firefly(target)

# For a list of targets
targets = ('WASP-43 b', 'WASP-18 b', 'WASP-91 b', 'WASP-12 b',
           'WASP-126 b', 'LHS 3844 b', 'GJ 1252 b', 'TOI-270 b')
firefly(targets)

# For a list of random viable tess targets with full prior info
targets = tess_viable(k=10)
firefly(targets)

# To simply use the default fitting settings on a random viable target
firefly()

```
Below shows a working example with all variables.
```python
from firefly import firefly, tess_viable

def main():
    targets = tess_viable(k=10)
    firefly(
        # Firefly Interface
        targets, 
        archive='nasa', 
        curve_sample=0.01, 
        email=True,
        to=['your.email@server.com'], 
        clean=False,
        cache=False,
        # TransitFit Variables
        cutoff=0.25,
        window=2.5,
        nlive=300, 
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
    )

main()
```
#### Console Output
```
Target search chose HD 2685 b.
Input checks passed.

Priors generated from the NASA Archive for HD 2685 b (TIC 267263253).

+-------------+----------------+----------------+------------+----------+
| Parameter   | Distribution   |        Input_A |    Input_B | Filter   |
|-------------+----------------+----------------+------------+----------|
| P           | gaussian       |    4.12689     |   0.004126 |          |
| t0          | gaussian       |    2.45868e+06 |   0.02345  |          |
| a           | gaussian       |    0.0568      |   0.0006   |          |
| inc         | gaussian       |   89.252       |   0.415    |          |
| rp          | uniform        |    0.0469772   |   0.187909 | 0        |
| host_T      | fixed          | 6844.5         | 140        |          |
| host_z      | fixed          |    0.02        |   0.06     |          |
| host_r      | fixed          |    1.575       |   0.055    |          |
| host_logg   | fixed          |    4.22        |   0.195    |          |
+-------------+----------------+----------------+------------+----------+

Splitting the lightcurve into seperate epochs using the following parameters.

+-------------+---------------+
| Parameter   |         Value |
|-------------+---------------|
| t0          |   2.45868e+06 |
| P           |   4.12689     |
| t14         | 265.74        |
+-------------+---------------+

Searching MAST for HD 2685 b.

Query from MAST returned 3 data products for HD 2685 b (TIC 267263253).

Downloading https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018206045859-s0001-0000000267263253-0120-s_lc.fits
|==========================================| 2.0M/2.0M (100.00%)         1s
Downloading https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2020186164531-s0027-0000000267263253-0189-s_lc.fits
|==========================================| 1.7M/1.7M (100.00%)         1s
Downloading https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2020212050318-s0028-0000000267263253-0190-s_lc.fits
|==========================================| 1.8M/1.8M (100.00%)         1s

+----------------+---------------------------------------------------------+
| Sector         | Product                                                 |
|----------------+---------------------------------------------------------|
| TESS Sector 1  | tess2018206045859-s0001-0000000267263253-0120-s_lc.fits |
| TESS Sector 27 | tess2020186164531-s0027-0000000267263253-0189-s_lc.fits |
| TESS Sector 28 | tess2020212050318-s0028-0000000267263253-0190-s_lc.fits |
+----------------+---------------------------------------------------------+

Splitting up the lightcurves into seperate epochs:

Light curve 6 discarded
Light curve 3 discarded
Light curve 6 discarded

A total of 18 lightcurves were generated.

A random sample of 18 lightcurves will be fitted across all TESS Sectors.
```

Parameters
----------
    
#### targets : str, list
A list of exoplanet targets.
Input is a list tuple of strings:
```python
('WASP-43 b', 'WASP-18 b', 'WASP-91 b')
```    

#### archive : {`'eu'`, `'nasa'`}, optional
The exoplanet archive to use for priors. Supports:

- 'eu'
- 'nasa'
The default is `'eu'`.

#### curve_sample : int {0 < curve_sample <= 1}, optional
The fraction of curves generated to fit against. For example, setting
```python
curve_sample = 0.5
```
will fit half the curves extracted. The formula for this works as:
```python
total_curves = curves_extracted * curve_sample
```
Always returns an int. For example:
```python
curve_sample = 0.001
```
will fit using only 1 lightcurve from each sector. 
The default is `1` to fit all lightcurves across all sectors.

#### email : bool, optional
If True will send status emails. The default is `False`.

#### to : str, optional
The email address to send status updates to.
```python
to=['transitfit.server@gmail.com']
```

#### clean : bool, optional
If True will delete all downloaded files and zip outputs only.
The default is `False`.

#### cutoff : float, optional
If there are no data within 
```python
t14 * cutoff of t0, 
```
a period will be discarded. Default is `0.25`.

#### window : float, optional
Data outside of the range 
```python
[t0 ± (0.5 * t14) * window] 
```
will be discarded.

#### nlive : int, optional
The number of live points to use in the nested sampling retrieval.
Default is `1000`.

#### detrending_list : array_like, shape (n_detrending_models, 2)
A list of different detrending models. Each entry should consist
of a method and a second parameter dependent on the method.
Accepted methods are:
```python
['nth order', order]
['custom', function, [global fit indices, filter fit indices, epoch fit indices]]
['off', ]
```
Function here is a custom detrending function. TransitFit assumes
that the first argument to this function is times and that all
other arguments are single-valued - TransitFit cannot fit
list/array variables. If `'off'` is used, no detrending will be
applied to the **LightCurves** using this model.
If a custom function is used, and some inputs to the function
should not be fitted individually for each light curve, but should
instead be shared either globally, within a given filter, or within
a given epoch, the indices of where these fall within the arguments
of the detrending function should be given as a list. If there are
no indices to be given, then use an empty list: []
e.g. if the detrending function is given by
```python
foo(times, a, b, c):
    do something
```
and a should be fitted globally, then the entry in the method_list
would be 
```python
['custom', foo, [1], [], []].
```

#### dynesty_sample : {`'unif'`, `'rwalk'`, `'rstagger'`, `'slice'`, `'rslice'`, `'hslice'`}, optional
Method used to sample uniformly within the likelihood constraint,
conditioned on the provided bounds. Unique methods available are:

- uniform sampling within the bounds(`'unif'`) 
- random walks with fixed proposals (`'rwalk'`) 
- random walks with variable ("staggering") proposals (`'rstagger'`) 
- multivariate slice sampling along preferred orientations (`'slice'`) 
- "random" slice sampling along all orientations (`'rslice'`) 
- "Hamiltonian" slices along random trajectories (`'hslice'`) 
and any callable function which follows
the pattern of the sample methods defined in dynesty.sampling.
'auto' selects the sampling method based on the dimensionality of
the problem (from ndim). 
- When ndim < 10, this defaults to `'unif'`.
- When 10 <= ndim <= 20, this defaults to `'rwalk'`. 
- When ndim > 20, this defaults to `'hslice'` if a gradient is provided 
  and `'slice'` otherwise. 
- `'rstagger'` and `'rslice'` are provided as alternatives for
  `'rwalk'` and `'slice'`, respectively. Default is `'auto'`.

#### fitting_mode : {`'auto'`, `'all'`, `'folded'`, `'batched'`}, optional
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
Default is `'auto'`.

#### fit_ttv : boolean, optional
DESCRIPTION. The default is `False`.

#### limb_darkening_model : {`'linear'`, `'quadratic'`, `'squareroot'`, `'power2'`, `'nonlinear'`}, optional
The limb darkening model to use. Allowed models are

- `'linear'`
- `'quadratic'`
- `'squareroot'`
- `'power2'`
- `'nonlinear'`
With the exception of the non-linear model, all models are constrained
by the method in Kipping (2013), which can be found at
https://arxiv.org/abs/1308.0009. Use `ldc_low_lim` and `ldc_high_lim`
to control the behaviour of unconstrained coefficients.
Default is `'quadratic'`.

#### ld_fit_method : {`'coupled'`, `'single'`, `'independent'`, `'off'`}, optional
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

#### max_batch_parameters : int, optional
The maximum number of parameters to use in a single retrieval.
Default is `25`.

#### batch_overlap : int, optional
The number of epochs to overlap in each batch. This will be adhered
to where possible. Default is `2`.

#### dlogz : float, optional
Retrieval iteration will stop when the estimated contribution of
the remaining prior volume to the total evidence falls below this
threshold. Explicitly, the stopping criterion is
```python
ln(z + z_est) - ln(z) < dlogz,
```
where z is the current evidence
from all saved samples and z_est is the estimated contribution from
the remaining volume. The default is 
```python
1e-3 * (nlive - 1) + 0.01.
```

#### maxiter : int or `None`, optional
The maximum number of iterations to run. If `None`, will
continue until stopping criterion is reached. Default is `None`.

#### maxcall : int or `None`, optional
The maximum number of likelihood calls in retrieval. If None, will
continue until stopping criterion is reached. Default is `None`.

#### dynesty_bounding : {`'none'`, `'single'`, `'multi'`, `'balls'`, `'cubes'`}, optional
The decomposition to use in sampling. Default is `'multi'`.

#### normalise : bool, optional
If True, will assume that the light curves have not been normalised and
will fit normalisation constants within the retrieval. The range to
fit normalisation constants c_n are automatically detected using
```python
1/f_min <= c_n <= 1/f_max
```
as the default range, where f_min and f_max are the minimum and maximum
flux values for a given light curve. Default is `True`.

#### detrend : bool, optional
If True, will initialise detrending fitting. Default is `True`.

Returns
-------
A whole lot of data to science!
Zipped files are found in:
```
firefly/WASP-43 b timestamp.gz.tar
```
