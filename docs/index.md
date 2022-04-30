## Firefly

# 🛰 **firefly**
### A target selector for use with TransitFit to fit TESS lightcurves.

![console](https://user-images.githubusercontent.com/10788239/147519786-e4e1e856-9dca-4350-947b-fc5c16b43763.gif)

<p align="center">
  <img src="https://github.com/sourestdeeds/firefly/blob/main/firefly/data/WASP-100%20b%20density.png">
</p>

#### Installation
```bash
pip install firefly-tess
```
or
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
```python
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
```python
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

priors('wasp190b', 'spearnet')
```
```python
Priors generated from the SPEARNET Archive for WASP-190 b (TIC 116156517).

+-------------+----------------+---------------+---------------------+----------+
| Parameter   | Distribution   |       Input A | Input B             | Filter   |
|-------------+----------------+---------------+---------------------+----------|
| P           | fixed          |    5.36773    |                     |          |
| t0          | gaussian       |    2.4578e+06 | 0.007               |          |
| a           | gaussian       |    0.0643173  | 0.0032913001497495  |          |
| inc         | gaussian       |   86.5471     | 0.1529545135471545  |          |
| w           | fixed          |   90          |                     |          |
| ecc         | fixed          |    0          |                     |          |
| rp          | uniform        |    0.0730344  | 0.08926423524371321 | 0        |
| t14         |                |  286.121      |                     |          |
| host_T      |                | 6400          | 128.0               |          |
| host_z      |                |   -0.02       | 0.05                |          |
| host_r      |                |    1.6        | 0.08000000000000002 |          |
| host_logg   |                |    4.17       | 0.0834              |          |
+-------------+----------------+---------------+---------------------+----------+
```

# Citing Firefly

```
@misc{firefly,
  author = {Stephen Charles and Joshua Hayes},
  title = {Firefly},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/sourestdeeds/firefly}},
}
```
