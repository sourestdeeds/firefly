#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 13:14:20 2021

@author: sourestdeeds
"""

import pandas as pd
import os
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns; sns.set_theme(style="whitegrid")


def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value
        patch.set_width(new_value)
        patch.set_x(patch.get_x() + diff * .5)


def plot_epoch():
    fig, ax = plt.subplots(figsize=(20,25))
    archive_list = ['eu', 'nasa', 'org', 'all']
    for i, archive in enumerate(archive_list):
        here = os.path.dirname(os.path.abspath(__file__))
        df = pd.read_csv(f'{here}/data/Targets/{archive}_tess_viable.csv') \
        .sort_values('Epochs', ascending=False).reset_index(drop=True)
        # df = pd.read_csv(f'{archive}_tess_viable.csv') \
        # .sort_values('Epochs', ascending=False).reset_index(drop=True)
        candidates = f'{len(df)} {archive.upper()} Archive ' +\
                     'Candidates Ranked by Epoch Count ' +\
                     f"({df['Epochs'].sum()} Total)"
        df[candidates] = 1
        df['Candidate Cumsum'] = df[candidates].cumsum()
        total_candidates = df[candidates].sum()
        df['Candidate Frequency'] = df['Candidate Cumsum'] / total_candidates
        
        df['Epoch Cumsum'] = df['Epochs'].cumsum()
        total_epochs = df['Epochs'].sum()
        df['Epoch Frequency'] = df['Epoch Cumsum'] / total_epochs
        # highlights = df.groupby('Frequency')['Exoplanet', 'Epochs']
        # .agg({'Exoplanet':'count','Epochs':'sum'})

        highlights = df.groupby(pd.cut(df["Epoch Frequency"], 
                                       np.arange(0, 1.0+0.1, 0.1))).sum()
        cumsum_cand = highlights[candidates].cumsum()
        temp = 411 + i
        ax=plt.subplot(temp)
        [ax.axvline(x=i+0.45, color='k', marker=',', 
                    alpha=0.5, ymin=0,ymax=0.25) for i in range(9)]
        [ax.axvline(x=i+0.45, color='k', marker=',', 
                    alpha=0.5, ymin=0.75,ymax=1) for i in range(9)]
        [ax.text(i+0.4,total_epochs/30, 
                 f'{i+1}0% - {cumsum_cand[i]} Candidates',rotation=90) 
                 for i in range(10)]
        sns.barplot(ax=ax, data=highlights, 
                    x=candidates, y = 'Epochs', 
                    hue='Epochs', dodge=False)
        change_width(ax, 0.7)
        plt.legend(title='Epochs',loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.show()
  
