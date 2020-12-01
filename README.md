# target

## Dependancies
        lightkurve
        transitfit



Generates the priors and host star variables for a chosen target.
Requires a sub-directory /data/ to store exoplanet csv files from the archive.
Target lightcurve files must be contained within a folder by the name of its target.
Example use with TransitFit is the following:
        
        # Host Info
        exoplanet, curves = 'WASP-43 b', 1
        host_T, host_z, host_r, host_logg  = target(exoplanet, curves)
        
        # Paths to data, priors, and filter info:
        data = 'data/data_paths.csv'
        priors = 'data/priors.csv'
        
    Parameters
    ----------
    exoplanet : string, example: 'WASP-43 b'
        The target exoplanet.
        
    curves : int, optional
        How many light curves to fit. Updates data paths for chosen target.
        Must be contained within the target folder, ie 'WASP-43 b/split_curve_0.csv'.
        The default is 1.
        
    dtype : string, optional
        Allows for inputs 'nasa' or 'eu'. The default is 'eu'.
        
        EU : File must be renamed and contained within directory 'data/eu_data.csv'.
        http://exoplanet.eu/catalog/#
        
        NASA : File must be renamed and contained within directory 'data/nasa_data.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html
        
    Returns
    -------
    host_T : tuple or None
        The effective temperature of the host star, in Kelvin. 
    host_z : tuple or None
        The log_10 of the surface gravity of the host star, with gravity measured in cm/s2. 
    host_r : tuple or None
        The metalicity of the host, given as a (value, uncertainty) pair.
    host_logg : tuple or None
        The host radius in Solar radii.
    data/data_paths.csv : file
        The locations of the light curves to fit.
    data/priors.csv : file
        The priors for the chosen target.
   
