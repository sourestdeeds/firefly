# 🛰 **Firefly**
### A target selector for use with TransitFit to fit TESS lightcurves.

The focus of this investigation is Transit Photometry, which was first suggested as a
method of exoplanet discovery by Struve (1952), and expanded upon by Rosenblatt
(1971). HD 209458 b was the first to be found via the transit method (Charbonneau
et al., 2000), which was previously confirmed via Radial Velocity measurements (Mazeh
et al., 2000). Combined detections such as these allow both the minimum planet mass
and radius to be measured, which can then be compared with planetary evolution
predictions (Guillot, 2005; Baraffe et al., 2008) to infer the majority of planetary
ephemera. For this reason, confirmed exoplanets detected via multiple detection
methods, provide a more complete and accurate constraints on the parameters.

The Transiting Exoplanet Survey Satellite (Ricker et al., 2014) (TESS) is an all-sky
transit survey, whose primary goal is to detect Earth-sized planets orbiting bright stars,
allowing follow-up observations to determine planet masses and atmospheric
compositions. TESS has an 85% sky coverage, of which each sector is continuously
observed for 4 weeks. For higher ecliptic lattitudes, the sectors overlap creating
photometric time series for durations up to a year. The upper and lower ecliptic poles
are called the continuous viewing zones (CVZ), and are constantly
observed in a yearly rotation between the two poles regardless of sector. Such
multi-sector photometry allows for a steady stream of transits to be observed, which
lends itself well to probe for transit timing variations (TTV’s). Increasing the accuracy
of known parameters through the use of lightcurve fitting programs also benefits from a
consistent single source of observations, as the systematic variance between sectors is
minimal. TESS aims for 50 ppm photometric precision on stars with a TESS magnitude
of 9-15.

In this work, I make use of a novel transit fitting program, TransitFit (Hayes et al.,
2021a), which is capable of using information about the host and planet parameters,
alongside the observation filters to couple stellar limb-darkening coefficients across
wavelengths. It was primarily designed for use with transmission spectroscopy studies,
and employs transit observations at various wavelengths from different telescopes to
simultaneously fit transit parameters using nested sampling retrieval.

![console](https://user-images.githubusercontent.com/10788239/147519786-e4e1e856-9dca-4350-947b-fc5c16b43763.gif)

<p align="center">
  <img src="https://raw.githubusercontent.com/sourestdeeds/firefly/main/firefly/data/WASP-100%20b%20density.webp#center">
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
# Mast

The Mikulski Archive for Space Telescopes (MAST) is a NASA-funded project to
support, and provide to the public, a variety of astronomical archives, with a primary
focus on space-based telescopes operating in the optical, ultraviolet, and near-infrared
parts of the spectrum.

The time series photometry (Tenenbaum and Jenkins, 2018) captured by TESS was
mainly at 120 second cadence for its primary mission, with baselines ranging from
approximately 27 days to a year depending on sector overlapping. The entirety of CCD’s
are known as Full-Frame Images (FFI). Groups of pixels are downloaded at shorter cadence
to obtain a faster cadence for a subset of targets, known as Target Pixel (TP) files. Pixels
around the star are stored as arrays in the TP, one image per stamp. Aperture photometry
is then performed on each image which creates an array of fluxes known as Light Curve
(LC) files.

Lightcurves contain flux time series derived from calibrated two minute and twenty
second target pixels. Two photometry types are available, Simple Aperture Photometry
(SAP) and Pre-Search Data Conditioning SAP (PDCSAP). PDCSAP is a flux time
series which has common instrumental systematics removed using the Cotrending Basis
Vectors (CBV). The CBV’s represent a set of systematic trends present in the lightcurve
data for each CCD.
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
# Priors

To establish a priority list two conditions must be considered. The minimum
ephemera required for retrieval consist of the period P, transit mid-point t0,
semi-major axis a, inclination i and the ratio of the planet to its host star rp. The
additional values of eccentricity ecc and periastron w are optional, and taken into
account if they exist. Otherwise default values of ecc 0 and w 90 are used.
A viable candidate should satisfy both the existence of a full set of ephemera and a light
curve data product. In lieu of this, each archive was tested to ascertain which provided
the most candidates with full prior information, comprising of the Extrasolar Planets
Encyclopedia (Schneider et al., 2011), the NASA Exoplanet Archive (Akeson et al.,
2013), and the Exoplanet Orbit Database (Wright et al., 2011).

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
  author = {Stephen Charles and Joshua Hayes and Eamonn Kerins},
  title = {Firefly - A target selector for use with TransitFit to fit TESS lightcurves.},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/sourestdeeds/firefly}},
}
```
