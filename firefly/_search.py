#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spelling checker for the search engine.

@author: Steven Charles-Mindoza
"""


from ._archive import _download_eu, _download_nasa
from ._utils import suppress_print

from fuzzywuzzy import process
from pandas import read_csv


def _fuzzy_search(exoplanet, archive='eu'):
    if archive == 'eu':
        _download_eu()
        eu_csv = 'firefly/data/eu.csv'
        exo_list = read_csv(eu_csv, usecols=['# name']).values.tolist()
        exo_list = [j for i in exo_list for j in i]
    elif archive == 'nasa':
        _download_nasa()
        nasa_csv = 'firefly/data/nasa.csv'
        exo_list = read_csv(nasa_csv, usecols=['pl_name']).values.tolist()
        exo_list = [j for i in exo_list for j in i]
    ratios = process.extract(exoplanet,exo_list)
    highest = process.extractOne(exoplanet,exo_list)
    return highest, ratios

