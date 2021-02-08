# **firefly**
A target selector for use with TransitFit to fit TESS lightcurves.

#### Installation
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
astroquery
natsorted
```

- Targets passed are corrected for basic user input. 'wasp43b' is
interpreted as 'WASP-43 b'. List must be of the form given in the example below.
- Initial checks for targets from the exoplanet archive are then taken to ascertain 
whether the prior data extracted has entries in all columns. If there are missing
entries for a given target in the list, the user will be asked whether to proceed.
- All available split curves are fitted with TransitFit, then the results
are then zipped up and time stamped. Optionally sends an email upon an error or 
full completion of a target.

Background tasks for feeding data to TransitFit include:
- Set the filter path to the TESS Filter.
- Download EU/NASA/ORG/OEC exoplanet archive data every 2 days (checks the file age).
- Download MAST lightcurve data for target TESS Sectors.
- Split the lightcurves into seperate transits or epochs.
- Create the data paths to each seperate epoch.
- Run TransitFit.
- Delete all downloaded data and inputs.
- Zip up and timestamp the output.

An example use with TransitFit is the following:
```python
from firefly import firefly

firefly('wasp43b')
```
Below shows a working example with all variables.
```python
from firefly import firefly, tess_viable

    firefly(
            # Firefly Interface
            exoplanet,
            archive='eu',
            curve_sample=1,
            clean=False,
            cache=False,
            auto=True,
            # MAST Search
            hlsp=['SPOC'],
            cadence=[120],
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
            detrend=True,
       )
```
#### Console Output
```
Target search chose WASP-43 b.

Priors generated from the EU Archive for WASP-43 b (TIC 36734222).

+-------------+----------------+----------------+---------------+----------+
| Parameter   | Distribution   |        Input A |       Input B | Filter   |
|-------------+----------------+----------------+---------------+----------|
| P           | gaussian       |    0.813478    |   8.13478e-05 |          |
| t0          | gaussian       |    2.45573e+06 |   0.007       |          |
| a           | gaussian       |    0.01526     |   0.0001526   |          |
| inc         | gaussian       |   82.33        |   0.8233      |          |
| w           | gaussian       |  328           | 164           |          |
| ecc         | gaussian       |    0.0035      |   0.00175     |          |
| rp          | uniform        |    0.143652    |   0.175575    | 0        |
| host_T      | fixed          | 4520           |  90.4         |          |
| host_z      | fixed          |   -0.01        |   0.012       |          |
| host_r      | fixed          |    0.667       |   0.03335     |          |
| host_logg   | fixed          |    4.64544     |   0.00971112  |          |
+-------------+----------------+----------------+---------------+----------+

Splitting the lightcurve into seperate epochs using the following parameters.

+-------------+--------------+
| Parameter   |        Value |
|-------------+--------------|
| t0          |  2.45573e+06 |
| P           |  0.813478    |
| t14         | 87.89        |
+-------------+--------------+

Searching MAST for WASP-43 b (TIC 36734222).

Query from MAST returned 1 data products for WASP-43 b (TIC 36734222).

+-------------------------------------------------+----------+-----------+--------+-----------+
| Product                                         |   TIC ID |   Cadence | HLSP   | Mission   |
|-------------------------------------------------+----------+-----------+--------+-----------|
| tess2019058134432-s0009-0000000036734222-0139-s | 36734222 |       120 | SPOC   | TESS      |
+-------------------------------------------------+----------+-----------+--------+-----------+

```

Parameters
----------
    
#### targets : str, list
A list of exoplanet targets.
Input is a single string:
```python
'WASP-43 b'
```    

#### archive : {`'eu'`, `'nasa'`}, optional
The exoplanet archive to use for priors. Supports:

- 'eu'
- 'nasa'
- 'org'
- 'all'
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

#### clean : bool, optional
If True will delete all downloaded files and zip outputs only.
The default is `False`.

Returns
-------
A whole lot of data to science!
Zipped files are found in:
```
firefly/WASP-43 b timestamp.gz.tar
```
