# target

## Dependancies
        lightkurve
        transitfit

A target data retriever for confirmed/candidate TESS exoplanets.
Generates the priors and host star variables for a chosen target.
Downloads exoplanet archives every 10 days and stores in /data.
Target lightcurve files are downloaded from MAST, fits file is
stored in Planet/exoplanet/exoplanet.fits.
    
    An example use with TransitFit is the following:
        
    Host Info
        exoplanet = 'WASP-43 b'
        
        host_T, host_z, host_r, host_logg  = target(exoplanet)
        
    Paths to data, priors, and filter info:
        data = 'data/data_paths.csv'
        
        priors = 'data/priors.csv'
        
    Outputs
        results_output_folder = 'Planet/'+exoplanet+'/output_parameters'
        
        fitted_lightcurve_folder = 'Planet/'+exoplanet+'/fitted_lightcurves'
    
        plot_folder = 'Planet/'+exoplanet+'/plots'
    Parameters
    ----------
    exoplanet : string, example: 'WASP-43 b'
        The target exoplanet.
        
    curves : int, optional
        How many light curves to fit. Updates data paths for chosen target.
        Must be contained within the target folder, 
        ie 'WASP-43 b/split_curve_0.csv'.
        The default is 1.
        
    dtype : string, optional
        Allows for inputs 'nasa' or 'eu'. The default is 'nasa'.
        
        EU : Data downloaded and stored in 'data/eu_data.csv'.
        http://exoplanet.eu/catalog/#
        
        NASA : Data downloaded and stored in 'data/nasa_data.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html
    Returns
    -------
    host_T : tuple or None
        The effective temperature of the host star, in Kelvin. 
    host_z : tuple or None
        The log_10 of the surface gravity of the host star, 
        with gravity measured in cm/s2. 
    host_r : tuple or None
        The metalicity of the host, given as a (value, uncertainty) pair.
    host_logg : tuple or None
        The host radius in Solar radii.
    data/data_paths.csv : file
        The locations of the light curves to fit.
    data/priors.csv : file
        The priors for the chosen target.
    data/eu.csv : file
        EU : Data downloaded and stored in 'data/eu_data.csv'.
        http://exoplanet.eu/catalog/#
    data/nasa.csv : file
        NASA : Data downloaded and stored in 'data/nasa_data.csv'.
        https://exoplanetarchive.ipac.caltech.edu/index.html
    Planet/exoplanet/exoplanet.fits : file
        MAST : Data downloaded and stored in 'Planet/exoplanet/exoplanet.fits'.
        https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
