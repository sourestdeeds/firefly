#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spelling checker for the search engine.

@author: Steven Charles-Mindoza
"""


from ._archive import _download_nasa

from fuzzywuzzy import process
from pandas import read_csv


def _fuzzy_search(exoplanet):
    _download_nasa()
    nasa_csv = 'firefly/data/nasa.csv.gz'
    exo_list = read_csv(nasa_csv, usecols=['pl_name', 'tic_id']) \
               .dropna() .drop_duplicates('pl_name') \
               .drop(['tic_id'], axis=1) .values .tolist()
    exo_list = [j for i in exo_list for j in i]
    
    ratios = process.extract(exoplanet,exo_list)
    highest = process.extractOne(exoplanet,exo_list)
    return highest, ratios

