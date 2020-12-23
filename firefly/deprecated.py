#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Deprecated parallel auto_retrieval

@author: Steven Charles-Mindoza
"""
# Deprecated
'''
from .query import query
from ._utils import suppress_print, _iterable_target
from ._archive import _check_nan
from ._search import _fuzzy_search
from datetime import datetime
from shutil import rmtree, make_archive
from multiprocessing import Pool
from functools import partial
import sys
import os


def auto_retrieval(targets, processes=len(os.sched_getaffinity(0)) // 4,
                   archive='eu', curve_sample=1, nlive=300, detrending_list=[['nth order', 2]],
                   dynesty_sample='rslice', fitting_mode='folded', fit_ttv=False,
                   limb_darkening_model='quadratic', ld_fit_method='independent',
                   max_batch_parameters=25, batch_overlap=2, dlogz=None, 
                   maxiter=None, maxcall=None, dynesty_bounding='multi', 
                   normalise=True, detrend=True, email=False, 
                   to=['transitfit.server@gmail.com'], printing=False):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Input Checks
    if len(targets) == 1:
        processes = 1
        printing = True
    if not (archive == 'eu' or archive == 'nasa'):
        sys.exit('Archive data options for dtype are: \'eu\' or \'nasa\'')
    if not (0 < curve_sample <= 1):
        sys.exit('The curve sample must be in the range 0 < curve_sample <= 1.')
    exoplanet_list = []
    for i, exoplanet in enumerate(targets):
        highest, ratios = _fuzzy_search(exoplanet, archive=archive)
        exoplanet = highest[0]
        print(f'\nTarget search chose {highest[0]} from the list: \n')
        query(exoplanet, archive=archive)
        with suppress_print():  
            nan = _check_nan(exoplanet, archive=archive)
        if nan == True:
            _check_nan(exoplanet, archive=archive, printing=True)
            verify = ''
            while (verify!="y" and verify!="n"):
                verify = input(f'\nWARNING: {exoplanet} has missing '
                                   'prior entries. Proceed ([y]/n)?\n')
            if verify == "n":
                sys.exit()
            elif verify == "y":
                pass
        # If checks for nans are passed, continue
        exoplanet_list.append([exoplanet])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Parallel Processing
    func = partial(_iterable_target, 
                   archive=archive, email=email, printing=printing, nlive=nlive,
                   detrending_list=detrending_list, ld_fit_method=ld_fit_method,
                   dynesty_sample=dynesty_sample, fitting_mode=fitting_mode,
                   fit_ttv=fit_ttv, limb_darkening_model=limb_darkening_model, 
                   max_batch_parameters=max_batch_parameters, 
                   batch_overlap=batch_overlap, dlogz=dlogz, 
                   maxiter=maxiter, maxcall=maxcall, curve_sample=curve_sample,
                   dynesty_bounding=dynesty_bounding, 
                   normalise=normalise, detrend=detrend, to=to)
    try:
        with Pool(processes=processes) as pool:
            pool.map(func, exoplanet_list, chunksize=1)
    except KeyboardInterrupt:
        try:
            exo_folder = f'firefly/{exoplanet}'
            now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
            keyboard = f'firefly/KeyboardInterrupt/{exoplanet} ' +\
                       f'{now} KeyboardInterrupt'
            make_archive(keyboard, format='gztar',
                     root_dir=f'{os.getcwd()}/firefly/',
                     base_dir=f'{exoplanet}')
            rmtree(exo_folder)
        except:
            pass
    except BaseException:
        try:
            exo_folder = f'firefly/{exoplanet}'
            now = datetime.now().strftime("%d-%b-%Y %H:%M:%S")
            exception = f'firefly/Exception/{exoplanet} ' +\
                        f'{now} Exception'
            make_archive(exception, format='gztar',
                     root_dir=f'{os.getcwd()}/firefly/',
                     base_dir=f'{exoplanet}')
            rmtree(exo_folder)
        except:
            pass

 
def _iterable_target(exoplanet_list, archive='eu', curve_sample=1, nlive=300,
                     detrending_list=[['nth order', 2]],
                     dynesty_sample='rslice', fitting_mode='folded', fit_ttv=False,
                     limb_darkening_model='quadratic', ld_fit_method='independent',
                     max_batch_parameters=25, batch_overlap=2, dlogz=None, 
                     maxiter=None, maxcall=None, dynesty_bounding='multi', 
                     normalise=True, detrend=True, email=False, to=['transitfit.server@gmail.com'],
                     printing=False):
    for i, exoplanet in enumerate(exoplanet_list):
        try:
            # Printing suppressed within scope
            if printing == False:
                with suppress_print():
                    success = _retrieval(exoplanet, archive=archive, nlive=nlive,
                                         detrending_list=detrending_list,
                                         dynesty_sample=dynesty_sample,
                                         fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                                         limb_darkening_model=limb_darkening_model, 
                                         ld_fit_method=ld_fit_method,
                                         max_batch_parameters=max_batch_parameters, 
                                         batch_overlap=batch_overlap, dlogz=dlogz, 
                                         maxiter=maxiter, maxcall=maxcall, 
                                         dynesty_bounding=dynesty_bounding, 
                                         normalise=normalise, detrend=detrend,
                                         curve_sample=curve_sample)
            elif printing == True:
                success = _retrieval(exoplanet, archive=archive, nlive=nlive,
                                     detrending_list=detrending_list,
                                     dynesty_sample=dynesty_sample,
                                     fitting_mode=fitting_mode, fit_ttv=fit_ttv,
                                     limb_darkening_model=limb_darkening_model, 
                                     ld_fit_method=ld_fit_method,
                                     max_batch_parameters=max_batch_parameters, 
                                     batch_overlap=batch_overlap, dlogz=dlogz, 
                                     maxiter=maxiter, maxcall=maxcall, 
                                     dynesty_bounding=dynesty_bounding, 
                                     normalise=normalise, detrend=detrend,
                                     curve_sample=curve_sample)
                print(f'\nData location: {success}\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.')
            if email == True:
                _email(f'Success: {exoplanet}',
                       f'Data location: {success} \n\n'
                       'A new target has been fully retrieved across ' +
                       'all available TESS Sectors.', to=to)
        except Exception:
            trace_back = format_exc()
            if email == True:
                _email(f'Exception: {exoplanet}', trace_back, to=to)
            else:
                print(trace_back)
            pass
'''
