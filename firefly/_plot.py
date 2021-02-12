#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 13:14:20 2021

@author: Steven Charles-Mindoza
"""

import pandas as pd
import os
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
mpl.rcParams['figure.dpi'] = 300

base_context = {

                "font.size": 18,
                "axes.labelsize": 18,
                "axes.titlesize": 14,
                "xtick.labelsize": 16,
                "ytick.labelsize": 16,
                "legend.fontsize": 16,

                "axes.linewidth": 1.25,
                "grid.linewidth": 1,
                "lines.linewidth": 1.5,
                "lines.markersize": 6,
                "patch.linewidth": 1,

                "xtick.major.width": 0,
                "ytick.major.width": 0,
                "xtick.minor.width": 0,
                "ytick.minor.width": 0,

                "xtick.major.size": 0,
                "ytick.major.size": 0,
                "xtick.minor.size": 0,
                "ytick.minor.size": 0,

                }

sns.set_theme(style="whitegrid", rc = base_context)

def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value
        patch.set_width(new_value)
        patch.set_x(patch.get_x() + diff * .5)


def plot_epoch(sub=False):
    if sub==True:
        fig, ax = plt.subplots(figsize=(20,25))
    archive_list = ['eu', 'nasa', 'org', 'all']
    for i, archive in enumerate(archive_list):
        here = os.path.dirname(os.path.abspath(__file__))
        df = pd.read_csv(f'{here}/data/Targets/{archive}_tess_viable.csv') \
        .sort_values('Epochs', ascending=False).reset_index(drop=True)
        # df = pd.read_csv(f'Targets/{archive}_tess_viable.csv') \
        # .sort_values('Epochs', ascending=False).reset_index(drop=True)
        # candidates = f'{len(df)} {archive.upper()} Archive ' +\
        #            'Candidates Ranked by Epoch Count ' +\
        #             f"({df['Epochs'].sum()} Total)"
        candidates = f'{len(df)} {archive.upper()} Archive ' +\
             'Candidates Ranked by Descending Observed Transits'
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
        cumsum_epochs = highlights['Epochs'].cumsum()
        if sub==False:
            fig, ax = plt.subplots(figsize=(15,45))
            sns.despine(fig=fig, ax=ax)
        temp = 411 + i
        ax=plt.subplot(temp)
        [ax.axvline(x=i+0.4, color='k', marker=',',
                    alpha=0.5, ymin=0,ymax=0.15) for i in range(10)]
        [ax.axvline(x=i+0.45, color='k', marker=',',
                    alpha=0.5, ymin=0.8,ymax=1) for i in range(10)]
        # [ax.text(i+0.34,total_epochs/47,
        #          f"{i+1}0% - {cumsum_cand[i]} Candidates - {cumsum_epochs[i]} Transits",rotation=90)
        #          for i in range(10)]
        [ax.text(i+0.34,total_epochs/47,
                 f"{i+1}0% of All Transits - {round(cumsum_cand[i]*100/total_candidates)}% Candidates",rotation=90)
                 for i in range(10)]
        [ax.text(i-0.1, highlights['Epochs'][i] + 5,
                 str(highlights[candidates][i]),rotation=0) for i in range(10)]
        # Top Planets
        twenty_perc = cumsum_cand[2]
        textstr = '\n'.join(df['Exoplanet'][0:twenty_perc])
        props = dict(boxstyle='round', facecolor='white', alpha=0.1)
        ax.text(1.075, 0.96, f'30% of \nAll Transits \n{twenty_perc} Candidates',
                transform=ax.transAxes, fontsize=16, weight='bold',
                verticalalignment='top', bbox=props, ha='center')
        ax.text(1.02, 0.85, textstr, transform=ax.transAxes, fontsize=16,
                verticalalignment='top', bbox=props)
        # PLOT!
        highlights['Observed Transits'] = highlights['Epochs']
        sns.barplot(ax=ax, data=highlights, 
                    x=candidates, y = 'Observed Transits',
                    dodge=False, palette='rocket')
        change_width(ax, 0.6)
        # ax.get_legend().remove()
        
        ax.xaxis.tick_bottom()
        column_labels = [f'{i}0 - {i+1}0%'.replace('00 - 10%', '0 - 10%') for i in range(10)]
        ax.set_xticklabels(column_labels, minor=False)
        if sub ==False:
            plt.savefig(f'{archive}_epoch_rank.jpg', bbox_inches='tight')
    if sub==True:
        plt.savefig('epoch_rank.jpg', bbox_inches='tight')
        

def lc_plot(file):
    from lightkurve import LightCurve
    df = pd.read_csv(file).dropna()
    time = df['Time'].values - 2457000
    flux = df['Flux'].values
    flux_err = df['Flux err'].values
    
    lc = LightCurve(time, flux, flux_err, time_format='btjd')
    lc.remove_outliers()
    
    # pg = lc.normalize(unit='ppm') \
    #       .to_periodogram(oversample_factor=1)
    pg = lc.to_periodogram(oversample_factor=1)
    condition = (pg.power==np.max(pg.power))
    #max_p = np.where(condition)[0]
    # pg = lc.normalize(unit='ppm') \
    #      .to_periodogram(oversample_factor=25, maximum_period=np.abs(max_p/60))
    pg = lc.to_periodogram(oversample_factor=25, maximum_period=np.abs(1.5))
    # PLOT!
    fig, ax = plt.subplots(figsize=(15,15))
    ax=plt.subplot(311)
    plt.errorbar(lc.time, lc.flux, 
                 lc.flux_err, color='b',
                 alpha=0.2, zorder=1, capsize=2, ls='none')
    plt.scatter(lc.time, lc.flux, color='k', s=0.5, alpha=0.6, zorder=2)
    plt.xlabel('Time (BTJD)')
    plt.ylabel('Flux')
    
    ax=plt.subplot(312)
    plt.plot(pg.period, pg.power, color='k', alpha=0.8)
    plt.xlabel('Period')
    plt.ylabel('Power')
    #pg.show_properties()

    # Folded
    period = pg.period_at_max_power
    condition = (lc.flux==np.min(lc.flux))
    t0 = np.where(condition)[0][0]
    t0 = lc.time[t0]
    folded_lc = lc.fold(period=period, t0=t0) 
    
    ax=plt.subplot(313)
    plt.errorbar(folded_lc.time, folded_lc.flux, 
                 folded_lc.flux_err, color='b',
                 alpha=0.2, zorder=1, capsize=2, ls='none')
    plt.scatter(folded_lc.time, folded_lc.flux, 
                alpha=0.4, color='k', zorder=2, s=2)
    plt.xlabel('Phase')
    plt.ylabel('Flux')
    filename = f'{file}'.replace('.csv', '')
    plt.savefig(f'{filename}.jpg')
