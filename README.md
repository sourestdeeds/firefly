# **firefly**
### A target selector for use with TransitFit to fit TESS lightcurves.
<p align="center">
  <img src="https://raw.githubusercontent.com/sourestdeeds/firefly/main/firefly/data/filter_0.png?token=ACSJ3D7C7KDFPAFUZD7RNULAK7E6A">
</p>

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

# firefly
```python
from firefly import firefly

firefly('wasp43b')
```
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
| w           | gaussian       |  328           |  20           |          |
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
# mast
```python
from firefly import mast

mast('wasp43b')
```
```
Searching MAST for WASP-43 b (TIC 36734222).

Query from MAST returned 3 data products for WASP-43 b (TIC 36734222).

+----------------------------------------------------------------+----------+-----------+----------+-----------+
| Product                                                        |   TIC ID |   Cadence | HLSP     | Mission   |
|----------------------------------------------------------------+----------+-----------+----------+-----------|
| hlsp_diamante_tess_lightcurve_tic-0000000036734222_tess_v1_llc | 36734222 |      1800 | DIAMANTE | TESS      |
| hlsp_qlp_tess_ffi_s0009-0000000036734222_tess_v01_llc          | 36734222 |      1800 | QLP      | TESS      |
| tess2019058134432-s0009-0000000036734222-0139-s                | 36734222 |       120 | SPOC     | TESS      |
+----------------------------------------------------------------+----------+-----------+----------+-----------+
```
# priors
```python
from firefly import priors

priors('wasp43b')
```
```
Priors generated from the EU Archive for WASP-43 b (TIC 36734222).

+-------------+----------------+----------------+---------------+----------+
| Parameter   | Distribution   |        Input A |       Input B | Filter   |
|-------------+----------------+----------------+---------------+----------|
| P           | gaussian       |    0.813478    |   8.13478e-05 |          |
| t0          | gaussian       |    2.45573e+06 |   0.007       |          |
| a           | gaussian       |    0.01526     |   0.0001526   |          |
| inc         | gaussian       |   82.33        |   0.8233      |          |
| w           | gaussian       |  328           |  20           |          |
| ecc         | gaussian       |    0.0035      |   0.00175     |          |
| rp          | uniform        |    0.143652    |   0.175575    | 0        |
| host_T      | fixed          | 4520           |  90.4         |          |
| host_z      | fixed          |   -0.01        |   0.012       |          |
| host_r      | fixed          |    0.667       |   0.03335     |          |
| host_logg   | fixed          |    4.64544     |   0.00971112  |          |
+-------------+----------------+----------------+---------------+----------+
```
