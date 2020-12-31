# firefly

A collection of tools to aid the user experience in the use of
TransitFit for confirmed TESS targets.

#### Installation (add later)
```bash
python setup.py install
```

#### Dependancies
```python
transitfit
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

```
