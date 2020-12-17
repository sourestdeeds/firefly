# firefly

A collection of tools to aid the user experience in the use of
TransitFit.

#### Installation (add later)
```bash
pip install firefly
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
target = 'WASP-43 b'
query(target)
```
```
SearchResult containing 1 data products.

 observation  target_name                     productFilename                    
------------- ----------- -------------------------------------------------------
TESS Sector 9    36734222 tess2019058134432-s0009-0000000036734222-0139-s_lc.fits

Priors generated from the EU Archive for WASP-43 b.

 Parameter Distribution       Input_A     Input_B Filter
         P     gaussian  8.134775e-01    0.070000       
        t0     gaussian  2.455727e+06    0.012000       
         a     gaussian  1.526000e-02    0.000180       
       inc     gaussian  8.233000e+01    0.200000       
        rp      uniform  7.980670e-02    0.319227      0
    host_T        fixed  4.520000e+03  120.000000       
    host_z        fixed -1.000000e-02    0.012000       
    host_r        fixed  6.670000e-01    0.010000       
 host_logg        fixed  4.645444e+00    0.019972       
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

Automated version of retrieval. Optionally sends an email upon an error or 
full completion of a target. Iteratively takes targets given and employs 
TransitFit across each TESS sector for every exoplanet in the list given.
If more than one target is given in a list, multiple cpu's will handle the extra
targets in seperate threads. The default for this behaviour is set to a
quarter of the maximum cores available to the current process.
All available split curves are fitted with TransitFit, then the results
are then zipped up and time stamped.

An example use with TransitFit is the following:
```python
from firefly import auto_retrieval

# For a single target
target = ('WASP-43 b',)
auto_retrieval(targets)

# For a list of targets
targets = ('WASP-43 b', 'WASP-18 b', 'WASP-91 b', 'WASP-12 b',
            'WASP-126 b', 'LHS 3844 b', 'GJ 1252 b', 'TOI-270 b')
auto_retrieval(targets)
```
