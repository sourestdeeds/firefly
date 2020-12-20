# firefly

A collection of tools to aid the user experience in the use of
TransitFit for confirmed TESS targets.

#### Installation (add later)
```bash
python setup.py install
```

#### Dependancies
```python
lightkurve
transitfit
```

## query

Performs a search from both MAST and the EU/NASA exoplanet
archive and prints the results to the console.

```python
from firefly import query
target = 'WASP-18 b'
query(target)
```
```
+----------------+--------------------------------------------------------------+----------+
| Observation    | Product                                                      |     Size |
|----------------+--------------------------------------------------------------+----------|
| TESS Sector 2  | tess2018234235059-s0002-0000000100100827-0121-s_lc.fits      |  2004480 |
| TESS Sector 3  | tess2018263035959-s0003-0000000100100827-0123-s_lc.fits      |  1998720 |
| TESS Sector 29 | tess2020238165205-s0029-0000000100100827-0193-s_lc.fits      |  1915200 |
| TESS Sector 30 | tess2020266004630-s0030-0000000100100827-0195-s_lc.fits      |  1998720 |
+----------------+--------------------------------------------------------------+----------+

Priors generated from the EU Archive for WASP-18 b.

+-------------+----------------+----------------+-------------+----------+
| Parameter   | Distribution   |        Input_A |     Input_B | Filter   |
|-------------+----------------+----------------+-------------+----------|
| P           | gaussian       |    0.941452    |   0.04      |          |
| t0          | gaussian       |    2.45422e+06 |   0.038     |          |
| a           | gaussian       |    0.02047     |   0.00038   |          |
| inc         | gaussian       |   86           |   2.5       |          |
| rp          | uniform        |    0.0486661   |   0.194664  | 0        |
| host_T      | fixed          | 6400           | 100         |          |
| host_z      | fixed          |    0           |   0.09      |          |
| host_r      | fixed          |    1.23        |   0.045     |          |
| host_logg   | fixed          |    4.35179     |   0.0347287 |          |
+-------------+----------------+----------------+-------------+----------+
 ```

## retrieval

A target data retriever for confirmed/candidate TESS exoplanets.
Generates the priors and host star variables for a chosen target.
Downloads exoplanet archives every 10 days and stores in /data.
Target lightcurve files are downloaded from MAST, then split into 
separate epochs. Upon user entry of the amount of epochs to fit,
TransitFit will fit the curves and return the results. The results
are then zipped up and time stamped.

An example use with TransitFit is the following:
```python
from firefly import retrieval
target = 'WASP-43 b'
retrieval(target)
```

## auto_retrieval

Automated version of retrieval.
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
from firefly import auto_retrieval

# For a single target
target = ('WASP-43 b',)
auto_retrieval(target)

# For a list of targets
targets = ('WASP-43 b', 'WASP-18 b', 'WASP-91 b', 'WASP-12 b',
            'WASP-126 b', 'LHS 3844 b', 'GJ 1252 b', 'TOI-270 b')
auto_retrieval(targets)
```
