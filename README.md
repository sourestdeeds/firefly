# firefly

A collection of tools to aid the user experience in the use of
TransitFit for confirmed TESS targets.

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
- If more than one target is given in a list, multiple cpu's will handle the extra
targets in seperate threads. The default for this behaviour is set to a
quarter of the maximum cores available to the current process.
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
