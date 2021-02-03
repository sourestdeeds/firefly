#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
                                        
                     (                   (     (        
                     )\ )  (   (      (  )\ )  )\ (     
                    (()/(  )\  )(    ))\(()/( ((_))\ )  
                     /(_))((_)(()\  /((_)/(_)) _ (()/(  
                    (_) _| (_) ((_)(_)) (_) _|| | )(_)) 
                     |  _| | || '_|/ -_) |  _|| || || | 
                     |_|   |_||_|  \___| |_|  |_| \_, | 
                                                  |__/  
                                  
                  A data retriever for use with TransitFit.

    Functions
    ---------
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    firefly()
    
        Automated version of retrieval. For a single target the procedure is:
         
         >>> from firefly import firefly
             target = ('WASP-43 b',)
             firefly(target)
         
        For a list of targets:
             
             >>> from firefly import firefly
                 targets = ('WASP-43 b', 'WASP-18 b')
                 firefly(targets)
         
        - Targets passed are corrected for basic user input. 'wasp43b' is
          interpreted as 'WASP-43 b'. List must be of the form given in the 
          example below.
        - Initial checks for targets from the exoplanet archive are then taken 
          to ascertain whether the prior data extracted has entries in all 
          columns. If there are missing entries for a given target in the list, 
          the user will be asked whether to proceed.
        - Iteratively takes the targets given and employs TransitFit across 
          each TESS sector for every exoplanet in the list given.
        - All available split curves are fitted with TransitFit, then the 
          results are then zipped up and time stamped. Optionally sends an 
          email upon an error or full completion of a target.
        - Email is also disabled by default. If enabled, status updates on 
          completion and exceptions with the full traceback are sent. 
          
        Background tasks for feeding data to TransitFit include:
        
        - Set the filter path to the TESS Filter.
        - Download EU/NASA exoplanet archive data every 10 days 
          (checks the file age).
        - Download MAST lightcurve data for target TESS Sectors.
        - Split the lightcurves into seperate transits or epochs.
        - Create the data paths to each seperate epoch.
        - Run TransitFit.
        - Delete all downloaded data and inputs.
        - Zip and timestamp the output.
            
        Input is capable of handling :
        
            >>> ('wasp43b', 'WASp18b', 'wasP91--b') etc
            
        Forces corrections based on classifier: 
        
            >>> 'WASP', 'LTT', 'GJ' etc
            
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
"""

name = 'firefly'
__version__ = '0.7.7'

from .auto_retrieval import firefly
from ._archive import priors, tess, gen_tess, gen_tess_viable_ttv
