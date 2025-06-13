#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 13:14:20 2021

@author: Stephen Charles
"""

import pandas as pd
import os
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import to_rgb
from matplotlib import colors
from astropy.timeseries import LombScargle
from astropy import units as u
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

sns.set_theme(style="whitegrid")
sns.set_context('paper')
sns.set_palette('deep')

def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value
        patch.set_width(new_value)
        patch.set_x(patch.get_x() + diff * .5)


def plot_epoch(sub=False):
    sns.set_theme(style="whitegrid", rc = base_context)
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
        

def lc_plot(file, flatten=False):
    sns.set_theme(style="whitegrid", rc = base_context)
    from lightkurve import LightCurve
    df = pd.read_csv(file).dropna()
    time = df['Time'].values - 2457000
    flux = df['Flux'].values
    flux_err = df['Flux err'].values
    
    lc = LightCurve(time, flux, flux_err, time_format='btjd')
    lc = lc.remove_outliers().normalize()
    if flatten==True:
        lc = lc.remove_outliers().flatten()
    
    pg = lc.to_periodogram(oversample_factor=25)
    period = float(pg.period_at_max_power/u.d)
    pg = lc.to_periodogram(oversample_factor=25, maximum_period=period*1.25,
                           minimum_period=period*0.75)
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
    plt.savefig(f'{filename}.jpg', bbox_inches='tight')


def mw():
    from astropy import units as u
    import astropy.coordinates as apycoords
    from mw_plot import MWSkyMap, MWPlot
    here = os.path.dirname(os.path.abspath(__file__))
    # Viable Targets
    df = pd.read_csv(f'{here}data/Targets/nasa_tess_viable.csv')
    ra = df['RA'].values * u.deg
    dec = df['DEC'].values * u.deg
    z = df['Epochs'].values
    distance = df['Distance'] .values * u.pc
    c = apycoords.SkyCoord(ra=ra, dec=dec, distance=distance, frame='icrs')
    
    # All Targets
    df2 = pd.read_csv('firefly/data/nasa.csv.gz')
    ra2 = df2['ra'].values * u.deg
    dec2 = df2['dec'].values * u.deg
    distance2 = df2['sy_dist'] .values * u.pc
    c2 = apycoords.SkyCoord(ra=ra2, dec=dec2, distance=distance2, frame='icrs')
    
    plot_instance = MWPlot(mode='face-on', center=(0, 0)*u.kpc, radius= 12*u.kpc,
                       unit=u.kpc, coord='galactic', annotation=True,  grayscale=False)
    plot_instance.fontsize = 35  # fontsize for matplotlib plotting
    plot_instance.figsize = (20, 20)  # figsize for matplotlib plotting
    plot_instance.dpi = 300  # dpi for matplotlib plotting
    plot_instance.cmap = 'hot'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
    #plot_instance.clim = (vmin, vmax) # colorbar range
    plot_instance.imalpha = 1  # alpha value for the milkyway image
    plot_instance.s = 100.0  # make the scatter points bigger
    plot_instance.tight_layout = True # whether plt.tight_layout() will be run
    

    plot_instance.mw_scatter(-c.galactic.cartesian.x,
                             c.galactic.cartesian.y, [z, 'Transits Observed'])
    # plot_instance.mw_scatter(-c2.galactic.cartesian.x,
    #                          c2.galactic.cartesian.y, [distance2, ''] )
    plot_instance.savefig('mw_zoom_out.png')
    
    
    plot_instance = MWPlot(mode='face-on', center=(0, 0)*u.kpc, radius= 2*u.kpc,
                       unit=u.kpc, coord='galactic', annotation=True,  grayscale=False)
    plot_instance.fontsize = 35  # fontsize for matplotlib plotting
    plot_instance.figsize = (20, 20)  # figsize for matplotlib plotting
    plot_instance.dpi = 300  # dpi for matplotlib plotting
    plot_instance.cmap = 'hot'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
    #plot_instance.clim = (vmin, vmax) # colorbar range
    plot_instance.imalpha = 1  # alpha value for the milkyway image
    plot_instance.s = 100.0  # make the scatter points bigger
    plot_instance.tight_layout = True # whether plt.tight_layout() will be run
    
    plot_instance.scatter(-c2.galactic.cartesian.x,
                             c2.galactic.cartesian.y, color='dimgrey', zorder=1)
    plot_instance.mw_scatter(-c.galactic.cartesian.x,
                             c.galactic.cartesian.y, [z, 'Transits Observed'],
                             zorder=2)
    plot_instance.savefig('mw.png')
    
    
    plot_instance = MWPlot(mode='face-on', center=(0, 0)*u.kpc, radius= 1*u.kpc,
                       unit=u.kpc, coord='galactic', annotation=True,  grayscale=False)
    plot_instance.fontsize = 35  # fontsize for matplotlib plotting
    plot_instance.figsize = (20, 20)  # figsize for matplotlib plotting
    plot_instance.dpi = 300  # dpi for matplotlib plotting
    plot_instance.cmap = 'hot'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
    #plot_instance.clim = (vmin, vmax) # colorbar range
    plot_instance.imalpha = 1  # alpha value for the milkyway image
    plot_instance.s = 100.0  # make the scatter points bigger
    plot_instance.tight_layout = True # whether plt.tight_layout() will be run
    
    plot_instance.scatter(-c2.galactic.cartesian.x,
                             c2.galactic.cartesian.y, color='dimgrey', zorder=1)
    plot_instance.mw_scatter(-c.galactic.cartesian.x,
                             c.galactic.cartesian.y, [z, 'Transits Observed'],
                             zorder=2)
    plot_instance.savefig('mw_zoom_in.png')
    

def oc(t0, t0_err, file='Complete_Results.csv', exoplanet=None):
    path_ttv = file
    data_ttv = pd.read_csv(path_ttv).set_index('Parameter') \
                .filter(like = 't0', axis=0).drop(['Telescope', 'Filter'], axis=1)
    
    t0_o = t0 # data['Best'] .values
    t0_oerr = t0_err # data['Error'].astype(float) .values
    t0_c = data_ttv['Best'] .values
    t0_cerr = data_ttv['Error'].astype(float) .values
    epoch_no = (data_ttv['Epoch'].astype(int) + 1) . values
    
    ominusc = t0_c - t0_o
    ominuscerr = t0_cerr - t0_oerr
    ominusc *= 24 * 60
    ominuscerr *= 24 * 60
    from sklearn.preprocessing import scale
    ominusc = scale(ominusc)

    
    # Generalised Lomb Scargle Floating Mean
    ls = LombScargle(epoch_no, ominusc, ominuscerr)
    frequency, power = ls.autopower(nyquist_factor=1, samples_per_peak=100)
    
    # False Alarm Probability
    fap = ls.false_alarm_probability(power.max())
    levels = [0.99, 0.5, 0.05, 0.01]
    false_alarm_levels = ls.false_alarm_level(levels)
    pack = {'FAP':false_alarm_levels, 'Percentage':levels}
    period = 1/frequency[np.argmax(power)]
    # Get the best frequency for plotting??
    fit_x = np.linspace(epoch_no.min(), epoch_no.max(), 1000)
    fit_y = ls.model(fit_x, frequency[np.argmax(power)])
    
    
    fig, ax = plt.subplots(figsize=(15,15))
    
    ax=plt.subplot(311)
        
    # Make figure and axes
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 1)

    oc_ax = fig.add_subplot(gs[0])
    ls_ax = fig.add_subplot(gs[1])

    #Plot data
    oc_ax.errorbar(epoch_no, ominusc, ominuscerr, marker='.',
                   elinewidth=0.8, color='dimgrey', linestyle='',
                   capsize=2, alpha=0.8, zorder=1)
    oc_ax.scatter(epoch_no, ominusc, marker='.', color='dimgrey', zorder=2)
    oc_ax.axhline(0, color='black', linestyle='--', linewidth=1)
    oc_ax.plot(fit_x, fit_y, color='red', alpha=0.8)
    
    #ls_ax.scatter(1/frequency, power, color='dimgrey', alpha=0.1, zorder=1, s=10)
    ls_ax.plot(1/frequency, power, linestyle='-', linewidth=0.75, zorder=2,
               marker='', color='k')
    pos = period * 1.97
    for (level, i) in zip(false_alarm_levels, levels):
        ls_ax.axhline(level, color='r', linestyle='--', alpha=0.8)
        ls_ax.annotate(f'{str(int(i*100))}'+'$\%$',
                       (pos, level*1.018), color='k', ha='center')
    ls_ax.annotate(f'Period: {int(period)}\n{fap*100:.2f}$\%$',
                        (period, power.max()*1.03), color='k', weight='bold', ha='center')
    # Sort out labels etc
    oc_ax.set_xlabel('Epoch')
    oc_ax.set_ylabel('O-C (minutes)')

    ls_ax.set_xlabel('Period (Epochs)')
    ls_ax.set_ylabel('Power')

    oc_ax.tick_params('both', which='both', direction='in', bottom=True, left=True)
    ls_ax.tick_params('both', which='both', direction='in', bottom=True, left=True)

    ls_ax.set_xscale('linear')
    ls_ax.set_xlim([0, period * 2])
    upper_y_fap = false_alarm_levels[2] * 1.5
    upper_y_pow = max(power) * 1.2
    upper_y = np.max([upper_y_fap, upper_y_pow])
    ls_ax.set_ylim([-0.01, upper_y])
    fig.tight_layout()
    if exoplanet==None:
        fig.savefig('O-C.jpg', bbox_inches='tight')
    else:
        fig.savefig(f"firefly/{exoplanet}/{exoplanet.lower().replace(' ', '')}_o-c.jpg",
                    bbox_inches='tight')

    
def make_epoch_converter(exoplanet):
    '''
    Converts between the epoch idx used in TF and the epoch idx used in
    file naming (they don't match because we removed some things)

    Returns dict with {epoch_idx : filename}
    '''
    data = pd.read_csv(f'firefly/{exoplanet}/data_paths.csv')
    #print(data.columns)
    converter = {}

    for i, row in data.iterrows():
        #print(row[' Epoch'])
        converter[row['Epochs']] = row['Path']#[15:]

    return converter

def find_epoch_no(epoch_idx, t0, P, epoch_converter):
    '''
    Given an epoch idx from TransitFit, finds the number of orbits that have
    passed to get to this epoch - Useful for things like O-C plots!
    '''
    # Load up the relevant data file
    path = epoch_converter[epoch_idx]

    data = pd.read_csv(path)

    # Get the last time value and calculate how many periods have elapsed between t0 and then.
    return (data['Time'].values[-1] - t0)// P
  
def oc_fold(t0, t0err, P, file='Complete_results.csv', exoplanet=None, longterm=True):
    '''
    Loads in the t0 and errors
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Read in TTV Data
    converter = make_epoch_converter(exoplanet)
    path_ttv = file
    data_ttv = pd.read_csv(path_ttv).set_index('Parameter') \
                .filter(like = 't0', axis=0).drop(['Telescope', 'Filter'], axis=1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Compute o-c
    t0_o = t0
    t0_oerr = t0err
    t0_c = data_ttv['Best'] .values
    t0_cerr = data_ttv['Error'].astype(float) .values
    epoch_no = np.array([find_epoch_no(i, t0_o, P, converter)
                         for i in range(len(t0_c))]) * P
    epoch_no = epoch_no - epoch_no[0]
    ominusc = t0_o  - t0_c
    ominuscerr = t0_cerr
    ominusc *= 24 * 60
    ominuscerr *= 24 * 60
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Chi 2
    from sklearn.metrics import mean_absolute_error
    loss = mean_absolute_error(ominusc,ominuscerr)
    dof = len(ominusc)
    chi2_red = np.sum( (ominusc)**2 / (t0_cerr)**2)/dof
    sigma = np.sqrt(2./len(ominusc))
    nsig = (chi2_red-1)/sigma
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Do the Lomb-Scargel stuff.
    ls = LombScargle(epoch_no, ominusc, ominuscerr)
                
    max_p = P * len(ominusc)//2
    if longterm==True:
        frequency, power = ls.autopower(nyquist_factor=1, samples_per_peak=75)
    else:
        frequency, power = ls.autopower(nyquist_factor=1, samples_per_peak=25,
                                    minimum_frequency=1/max_p)
    best_f = frequency[np.argmax(power)]
    best_P = 1/frequency[np.argmax(power)]
    epoch_phase = (epoch_no - (epoch_no //best_P) * best_P)/best_P
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # FAP
    fap = ls.false_alarm_probability(power.max())
    levels = [0.5, 0.05, 0.01]
    false_alarm_levels = ls.false_alarm_level(levels)
    pack = {'FAP':false_alarm_levels, 'Percentage':levels}
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Fold and Fits
    epoch_phase = (epoch_no  - (epoch_no //best_P) * best_P)/best_P
    fit_x = np.linspace(epoch_no.min(), epoch_no.max(), 10000)
    fit_y = ls.model(fit_x, frequency[np.argmax(power)])
    fit_x_phase = (fit_x - (fit_x //best_P) * best_P)/best_P
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Begin Plots
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(3, 1)
    plt.set_cmap('plasma')
    
    oc_ax = fig.add_subplot(gs[0])
    phase_ax = fig.add_subplot(gs[1])
    ls_ax = fig.add_subplot(gs[2])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # o-c Plot
    oc_ax.errorbar(epoch_no, ominusc, ominuscerr, marker='.',
                   elinewidth=0.8, color='dimgrey', linestyle='',
                   capsize=2, alpha=0.8, zorder=1)
    oc_ax.scatter(epoch_no, ominusc, marker='.', zorder=2,
                  c=epoch_no)
    oc_ax.axhline(0, color='black', linestyle='--', linewidth=1)
    oc_ax.plot(fit_x, fit_y, color='red', alpha=0.8)
    red = 'dof'
    from matplotlib.offsetbox import AnchoredText
    txt = AnchoredText(f'$\chi^2_{{{red}}} = {chi2_red:.2f}\, ({nsig:.2f}\sigma)$' +\
                        f'\n$\mu_{{error}}= {loss:.2f}$',
                        loc='upper right', frameon=False,
                        prop=dict(fontweight="bold"))
    oc_ax.add_artist(txt)
    oc_ax.set_ylim([np.nanmin(ominusc)*1.5, np.nanmax(ominusc)*2])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Phase Plot
    phase_ax.errorbar(epoch_phase, ominusc, ominuscerr, marker='.',
                      elinewidth=0.8, color='dimgrey', linestyle='', capsize=2,
                      alpha=0.8, zorder=1)
    phase_ax.scatter(epoch_phase, ominusc, c=epoch_no,
                     marker='.',  zorder=2, alpha=0.5)
    phase_ax.axhline(0, color='black', linestyle='--', linewidth=1)
    phase_ax.plot(fit_x_phase[np.argsort(fit_x_phase)],
               fit_y[np.argsort(fit_x_phase)], color='red', alpha=0.8)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # LS Plot
    ls_ax.plot(P/frequency, power, linestyle='-', linewidth=0.75, zorder=2,
               marker='', color='k')
    
    pos = P * best_P * 1.97 #len(epoch_no)*0.985
    for (level, i) in zip(false_alarm_levels, levels):
            ls_ax.axhline(level, color='r', linestyle='--', alpha=0.8)
            ls_ax.annotate(f'{str(int(i*100))}'+'$\%$',
                           (pos, level*1.018), color='k', ha='center')
    ls_ax.annotate(f'Peak: {P*best_P:.2f}\n{fap*100:.2f}$\%$',
                        (P*best_P, power.max()*1.03), color='k', weight='bold', ha='center')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Labels
    oc_ax.set_xlabel('Period (Days)')
    oc_ax.set_ylabel('O-C (Minutes)')
    phase_ax.set_xlabel('Phase')
    phase_ax.set_ylabel('O-C (Minutes)')

    ls_ax.set_xlabel('Period (Days)')
    ls_ax.set_ylabel('Power')

    oc_ax.tick_params('both', which='both', direction='in', bottom=True, left=True)
    phase_ax.tick_params('both', which='both', direction='in', bottom=True, left=True)
    ls_ax.tick_params('both', which='both', direction='in', bottom=True, left=True)
    
    ls_ax.set_xlim([0, P * best_P * 2])
    
    upper_y_fap = false_alarm_levels[2] * 1.5
    upper_y_pow = max(power) * 1.2
    upper_y = np.nanmax([upper_y_fap, upper_y_pow])
    ls_ax.set_ylim([-0.01, upper_y])
    ls_ax.set_xscale('linear')

    fig.tight_layout()
    if exoplanet==None:
        fig.savefig('O-C_fold.jpg', bbox_inches='tight')
    else:
        if longterm==True:
            fig.savefig(f"firefly/{exoplanet}/{exoplanet} o-c longterm.jpg",
                    bbox_inches='tight')
        else:
            fig.savefig(f"firefly/{exoplanet}/{exoplanet} o-c.jpg",
                    bbox_inches='tight')
    return chi2_red, nsig, loss, fap, best_P
    
  
def read_fitted_lc(exoplanet, transits):
    '''
    Reads in lc's and cleans the data.

    '''
    from transitfit.lightcurve import LightCurve
    epoch_no = transits
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Read in lc's
    df = pd.DataFrame(columns=(['Phase', 'Normalised flux',
                                'Flux uncertainty','Best fit curve']))
    for i in range(0,epoch_no):
        file = f'firefly/{exoplanet}/fitted_lightcurves/t0_f0_e{i}_detrended.csv'
        lc = pd.read_csv(file)[['Phase', 'Normalised flux',
                                'Flux uncertainty','Best fit curve']]
        df = pd.concat([df, lc], ignore_index=True)
    
    df = df.sort_values('Phase').reset_index(drop=True)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Zoom in on transits
    fit_combined = df['Best fit curve']
    transit_mask = np.ma.masked_values(fit_combined, 1.).mask
    transit_loc = np.where(transit_mask==False)
    window = int(len(transit_loc[0]) * 2.5)
    # Check if the transit is still in baseline of data, and there is a transit
    if window==0:
        start, stop = 0, len(fit_combined)
    try:
        if ( (transit_loc[0][0] + window > len(fit_combined)) 
            or (transit_loc[0][0] - window < 0) 
            or (transit_mask).all() ):
            start, stop = 0, len(fit_combined)
    except:
        pass
    else:
        start, stop = transit_loc[0][0] - window, transit_loc[0][-1] + window
        if start < 0:
            start, stop = 0, len(fit_combined)
    transit_mask[start:stop] = False
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Apply Mask
    df_transit_mask = df[~transit_mask]
    # Check if the data masked leaves a significant transit
    if len(df_transit_mask) < int(0.05*len(fit_combined)):
        df_transit_mask = df
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Sigma clipping
    from astropy.stats.sigma_clipping import sigma_clip
    from astropy.stats.funcs import mad_std
    
    fluxdiff = np.max(df['Best fit curve']) - np.min(df['Best fit curve'])
    if fluxdiff < 0.001:
        sigma = 4
    elif 0.001 < fluxdiff < 0.005:
        sigma = 5
    else:
        sigma = 6
    df_transit_mask['Residual'] = df_transit_mask['Normalised flux'] \
                                  - df_transit_mask['Best fit curve']
    clip = sigma_clip(df_transit_mask['Residual'], sigma=sigma, 
                      stdfunc=mad_std).mask
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Apply Mask
    df_clipped_mask = df_transit_mask[~clip]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Assign numpy arrays
    time_all = df_clipped_mask['Phase'].values
    flux_all = df_clipped_mask['Normalised flux'].values
    flux_err_all = df_clipped_mask['Flux uncertainty'].values
    fit_xall = df_clipped_mask['Phase'].values
    fit_yall = df_clipped_mask['Best fit curve'].values
    return time_all, flux_all, flux_err_all, fit_xall, fit_yall

def density_scatter(exoplanet, transits, P, cadence):
    """
    Scatter plot colored by 2d histogram
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.stats.funcs import mad_std
    from astropy.stats.sigma_clipping import sigma_clip
    from scipy.interpolate import interpn
    from scipy import stats
    time_all, flux_all, flux_err_all, fit_xall, fit_yall = read_fitted_lc(exoplanet, transits)
    x, y, yerr = time_all, flux_all, flux_err_all
    diff = flux_all - fit_yall
    mad = np.std(diff)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    from transitfit.lightcurve import LightCurve
    
    cadence_min = cadence / 60
    cad_bin = cadence_min / (P * 60 * 24)
    lc = LightCurve(x, y, yerr)
    obs_length = x.max() - x.min()
    n_bins = int((obs_length)/cad_bin)
    bins=[n_bins,n_bins]
    binned_phase, binned_flux, binned_err, binned_residuals = lc.bin(cad_bin, diff)
    madbin = np.std(binned_residuals)
    obs_depth = 1 - np.min(fit_yall)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1]) ),
                data, np.vstack([x,y]).T, method="splinef2d",
                bounds_error=False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z, diff = x[idx], y[idx], z[idx], diff[idx]
    cadences = len(y)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Histogram residuals
    unbinned_counts, bins = np.histogram(diff, bins=30)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Binned Histogram residuals
    binned_counts, _ = np.histogram(binned_residuals, bins)
    weighted_binned_counts = binned_counts * unbinned_counts.max()/binned_counts.max()
    from matplotlib.offsetbox import AnchoredText
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # All errorbars
    sns.set_theme(style="whitegrid")
    sns.set_context('paper')
    sns.set_palette('deep')
    txt = AnchoredText(f'$\sigma_{{unbinned}}={mad:.6f}$\n$' +\
                       f'\sigma_{{binned}}={madbin:.6f}$',
                       frameon=False,loc='lower right',
                       prop=dict(fontweight="bold"))
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2, height_ratios=[3,1], width_ratios=[10,1])
    ax = fig.add_subplot(gs[0])
    cb_ax = fig.add_subplot(gs[1])
    res_ax = fig.add_subplot(gs[2])
    hist_ax = fig.add_subplot(gs[3], sharey=res_ax)
    
    plt.set_cmap('hot')
    # Residuals
    res_ax.errorbar(x, diff, yerr, color='dimgrey',
                  alpha=0.1, zorder=1, capsize=2, ls='none')
    res_ax.scatter(x, diff, s=5, c=z, edgecolor='none', zorder=2)
    
    res_ax.scatter(binned_phase, binned_residuals,
                   s=5, alpha=1, color='white', zorder=4, edgecolor='none')
    res_ax.axhline(y=0, linestyle='--', color='k', zorder=5, lw=0.75)
    
    plt.set_cmap('hot')
    # Light Curve
    ax.errorbar(x, y, yerr, color='dimgrey',
                  alpha=0.1, zorder=1, capsize=2, ls='none')
    ax.scatter(x, y, c=z, zorder=2, s=5, edgecolor='none')
    ax.add_artist(txt)
    
    ax.scatter(binned_phase, binned_flux, zorder=4, s=5,
                color='white', edgecolor='none', alpha = 0.9)
    ax.plot(fit_xall, fit_yall, marker='', color='k', zorder=5, lw=1)
    # Histogram
    rgba_color = colors.to_rgba('dimgrey')
    facecolor = (rgba_color[0], rgba_color[1], rgba_color[2], 0.6)
    rgba_color = colors.to_hex(rgba_color, keep_alpha=True)
    hist_ax.hist(bins[:-1], bins, weights=unbinned_counts,
                 orientation='horizontal', color=facecolor,
                 edgecolor=rgba_color, histtype='stepfilled')
    hist_ax.axhline(y=0, linestyle='--', color='k', zorder=3, lw=0.75)
    
    rgba_color = colors.to_rgba('red')
    facecolor = (rgba_color[0], rgba_color[1], rgba_color[2], 0.6)
    rgba_color = colors.to_hex(rgba_color, keep_alpha=True)
    hist_ax.hist(bins[:-1], bins, weights=weighted_binned_counts,
                 orientation='horizontal', color=facecolor,
                 edgecolor=rgba_color, histtype='stepfilled',alpha=0.5)
    # Labels
    ax.xaxis.set_ticklabels([])
    cb_ax.axis('off')
    position = cb_ax.inset_axes([0.1,0.1,0.2,0.8])
    hist_ax.tick_params('both', which='both', direction='in',
                                 labelleft=False, labelbottom=False,
                                 right=False, top=False)
    res_ax.set_xlabel('Phase')
    ax.set_ylabel('Normalised Flux')
    res_ax.set_ylabel('Residual')
    from matplotlib.cm import ScalarMappable
    m = ScalarMappable(cmap='hot')
    cbar = plt.colorbar(m, ax=cb_ax, shrink=0.3, cax=position)
    cbar.set_label('Normalised Flux Density')
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    fig.savefig(f'firefly/{exoplanet}/{exoplanet} density.png',
                bbox_inches='tight')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # No errorbar
    sns.set_theme(style="darkgrid")
    sns.set_context('paper')
    sns.set_palette('deep')
    txt = AnchoredText(f'$\sigma_{{unbinned}}={mad:.6f}$\n$' +\
                       f'\sigma_{{binned}}={madbin:.6f}$',
                       frameon=False,loc='lower right',
                       prop=dict(fontweight="bold"))
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2, height_ratios=[3,1], width_ratios=[10,1])
    ax = fig.add_subplot(gs[0])
    cb_ax = fig.add_subplot(gs[1])
    res_ax = fig.add_subplot(gs[2])
    hist_ax = fig.add_subplot(gs[3], sharey=res_ax)
    
    plt.set_cmap('hot')
    res_ax.scatter(x, diff, s=5, alpha=0.8, c=z, edgecolor='none', zorder=1)
    res_ax.scatter(binned_phase, binned_residuals,
                   s=5, alpha=0.8, edgecolor='none', color='white', zorder=2)
    res_ax.axhline(y=0, color='k', linestyle='--', zorder=3, lw=0.75)
    # Lightcurve
    plt.set_cmap('hot')
    ax.xaxis.set_ticklabels([])
    ax.scatter(x, y, c=z, zorder=2, s=5, edgecolor='none')
    ax.plot(fit_xall, fit_yall, marker='', color='k', zorder=4, lw=1)
    ax.scatter(binned_phase, binned_flux, zorder=3, s=5, color='white',
               edgecolor='none')
    ax.add_artist(txt)
    # Histogram
    rgba_color = colors.to_rgba('dimgrey')
    facecolor = (rgba_color[0], rgba_color[1], rgba_color[2], 0.6)
    rgba_color = colors.to_hex(rgba_color, keep_alpha=True)
    hist_ax.hist(bins[:-1], bins, weights=unbinned_counts,
                 orientation='horizontal', color=facecolor,
                 edgecolor=rgba_color, histtype='stepfilled')
    hist_ax.axhline(y=0, linestyle='--', color='k', zorder=3, lw=0.75)
    
    rgba_color = colors.to_rgba('red')
    facecolor = (rgba_color[0], rgba_color[1], rgba_color[2], 0.6)
    rgba_color = colors.to_hex(rgba_color, keep_alpha=True)
    hist_ax.hist(bins[:-1], bins, weights=weighted_binned_counts,
                 orientation='horizontal', color=facecolor,
                 edgecolor=rgba_color, histtype='stepfilled',alpha=0.5)
   # Labels
    ax.xaxis.set_ticklabels([])
    cb_ax.axis('off')
    position = cb_ax.inset_axes([0.1,0.1,0.2,0.8])
    hist_ax.tick_params('both', which='both', direction='in',
                                 labelleft=False, labelbottom=False,
                                 right=False, top=False)
    res_ax.set_xlabel('Phase')
    ax.set_ylabel('Normalised Flux')
    res_ax.set_ylabel('Residual')
    from matplotlib.cm import ScalarMappable
    m = ScalarMappable(cmap='hot')
    cbar = plt.colorbar(m, ax=cb_ax, shrink=0.3, cax=position)
    cbar.set_label('Normalised Flux Density')
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    fig.savefig(f'firefly/{exoplanet}/{exoplanet} density noerr.png',
                bbox_inches='tight')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # No resid or err
    sns.set_theme(style="darkgrid")
    sns.set_context('paper')
    sns.set_palette('deep')
    txt = AnchoredText(f'$\sigma_{{unbinned}}={mad:.6f}$\n$' +\
                       f'\sigma_{{binned}}={madbin:.6f}$',
                       frameon=False,loc='lower right',
                       prop=dict(fontweight="bold"))
    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    
    plt.set_cmap('hot')
    ax.scatter(x, y, c=z, zorder=2, s=5, edgecolor='none')
    ax.scatter(binned_phase, binned_flux, zorder=3, s=5, color='white',
               edgecolor='none')
    ax.plot(fit_xall, fit_yall, marker='', color='k', zorder=3, lw=1)
    ax.add_artist(txt)
    plt.xlabel('Phase')
    ax.set_ylabel('Normalised Flux')
    fig.savefig(f'firefly/{exoplanet}/{exoplanet} density noresiderr.png',
                bbox_inches='tight')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # No resid
    sns.set_theme(style="whitegrid")
    sns.set_context('paper')
    sns.set_palette('deep')
    txt = AnchoredText(f'$\sigma_{{unbinned}}={mad:.6f}$\n$' +\
                       f'\sigma_{{binned}}={madbin:.6f}$',
                       frameon=False,loc='lower right',
                       prop=dict(fontweight="bold"))
    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    
    plt.set_cmap('hot')
    ax.errorbar(x, y, yerr, color='dimgrey',
                  alpha=0.1, zorder=1, capsize=2, ls='none')
    ax.scatter(x, y, c=z, zorder=2, s=5, edgecolor='none')
    ax.scatter(binned_phase, binned_flux, zorder=4, s=5, color='white',
               edgecolor='none')
    ax.plot(fit_xall, fit_yall, marker='', color='k', zorder=3, lw=1)
    # ax.errorbar(stat[1][:len(stat[0])], stat[0], stat_err[0]/6, color='yellow',
    #               alpha=0.7, capsize=2, ls='none', zorder=3)
    ax.add_artist(txt)
    plt.xlabel('Phase')
    ax.set_ylabel('Normalised Flux')
    fig.savefig(f'firefly/{exoplanet}/{exoplanet} density noresid.png',
                bbox_inches='tight')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    return mad, madbin, obs_depth, cadences


def plot_coupled_comparison():
    here = os.path.dirname(os.path.abspath(__file__))
    here = ''
    coupled = pd.read_csv(f'{here}data/spear.csv')
    uncoupled = pd.read_csv(f'{here}data/spear_uncoupled.csv').set_index('pl_name').drop('WD 1856+534 b', axis=0)
    uncoupled = uncoupled.reset_index()
    
    pl_name_coupled = coupled['pl_name']
    pl_name_uncoupled = uncoupled['pl_name']
    
    missing_uncoup = np.setdiff1d(pl_name_uncoupled, pl_name_coupled).tolist()
    missing_coup = np.setdiff1d(pl_name_coupled, pl_name_uncoupled).tolist()
    uncoupled = pd.read_csv(f'{here}data/spear_uncoupled.csv').set_index('pl_name').drop(missing_uncoup, axis=0).drop('WD 1856+534 b', axis=0)
    uncoupled = uncoupled.reset_index()
    coupled = pd.read_csv(f'{here}data/spear.csv').set_index('pl_name').drop(missing_coup, axis=0)
    coupled = coupled.reset_index()
    
    
    pl_rad_coupled = coupled['pl_orbsmax']
    pl_rad_uncoupled = uncoupled['pl_orbsmax']
    
    resid = pl_rad_coupled - pl_rad_uncoupled
    
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2, height_ratios=[2,2], width_ratios=[10,1])
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[2], sharex=ax)
    
    
    ax.plot(pl_rad_coupled, pl_rad_uncoupled, marker='.', lw=0)
    ax.plot(pl_rad_uncoupled, color='r', marker='', lw=0.5)    
    ax.set_xlabel('Exoplanet')
    ax.set_ylabel(r'$r_p/r_*$')
    ax.legend(['Coupled', 'Uncoupled'])
    
    ax2.plot(resid, marker='.', lw=0)
    ax2.axhline(y=0, linestyle='--', color='k', zorder=3, lw=0.75)
    ax2.set_xlabel('Exoplanet')
    ax2.set_ylabel(r'$r_p/r_*$ coupled - uncoupled')
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    
    
    return

def plot_ld_params():
    here = ''
    
    uncoupled_rank_sens = pd.read_csv(f'{here}data/spear_uncoupled.csv').sort_values(by='Sensitivity', ascending=False).set_index('pl_name').drop('WD 1856+534 b', axis=0)
    top10 = uncoupled_rank_sens[0:40].reset_index()['pl_name']
    
    coupled = pd.read_csv(f'{here}spear_coupled_ld_params.csv')
    uncoupled = pd.read_csv(f'{here}spear_uncoupled_ld_params.csv').set_index('pl_name').drop('WD 1856+534 b', axis=0)
    uncoupled = uncoupled.reset_index()
    
    pl_name_coupled = coupled['pl_name']
    pl_name_uncoupled = uncoupled['pl_name']
    
    missing_uncoup = np.setdiff1d(pl_name_uncoupled, pl_name_coupled).tolist()
    missing_coup = np.setdiff1d(pl_name_coupled, pl_name_uncoupled).tolist()
    uncoupled = pd.read_csv(f'{here}spear_uncoupled_ld_params.csv').set_index('pl_name').drop(missing_uncoup, axis=0).drop('WD 1856+534 b', axis=0)
    uncoupled = uncoupled.reset_index()
    coupled = pd.read_csv(f'{here}spear_coupled_ld_params.csv').set_index('pl_name').drop(missing_coup, axis=0)
    coupled = coupled.reset_index()
    coupled_filter = coupled['pl_name'].isin(top10)
    
    top10_coup = coupled[coupled_filter]
    top10_uncoup = uncoupled[coupled_filter]
    
    # LD PARAMS
    q0_coupled = coupled['q0']
    q0err_coupled = coupled['q0err']
    q1_coupled = coupled['q1']
    q1err_coupled = coupled['q1err']
    q0_uncoupled = uncoupled['q0']
    q0err_uncoupled = uncoupled['q0err']
    q1_uncoupled = uncoupled['q1']
    q1err_uncoupled = uncoupled['q1err']
    
    q0_coupled_top10 = top10_coup['q0']
    q0err_coupled_top10 = top10_coup['q0err']
    q1_coupled_top10 = top10_coup['q1']
    q1err_coupled_top10 = top10_coup['q1err']
    q0_uncoupled_top10 = top10_uncoup['q0']
    q0err_uncoupled_top10 = top10_uncoup['q0err']
    q1_uncoupled_top10 = top10_uncoup['q1']
    q1err_uncoupled_top10 = top10_uncoup['q1err']
    
    # PLOT
    fig = plt.figure(figsize=(16,8))
    gs = gridspec.GridSpec(4, 2, height_ratios=[2,2,1,1])
    
    # q0
    q0 = fig.add_subplot(gs[0])
    q0.plot(q0_coupled, marker='.', lw=0)
    q0.plot(q0_uncoupled, color='r', marker='.', lw=0)  
    
    q0.set_xlabel('Exoplanet')
    q0.set_ylabel(r'$q_0$')
    q0.legend(['Coupled', 'Uncoupled'])
    
    # q1
    q1 = fig.add_subplot(gs[2], sharex=q0)
    q1.plot(q1_coupled, marker='.', lw=0)
    q1.plot(q1_uncoupled, color='r', marker='.', lw=0) 
    
    q1.set_xlabel('Exoplanet')
    q1.set_ylabel(r'$q_1$')
    #q1.legend(['Coupled', 'Uncoupled'])
    
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    
    # resid
    residq0 = q0_coupled - q0_uncoupled
    residq1 = q1_coupled - q1_uncoupled
    resq0 = fig.add_subplot(gs[4], sharex=q0)
    resq0.plot(residq0, marker='.', lw=0, color='k')
    resq0.axhline(y=0, linestyle='--', color='k', zorder=3, lw=0.75)
    resq0.set_ylabel(r'$q_0$ c - u')
    
    resq1 = fig.add_subplot(gs[6], sharex=q0)
    resq1.plot(residq1, marker='.', lw=0, color='k')
    resq1.axhline(y=0, linestyle='--', color='k', zorder=3, lw=0.75)
    resq1.set_xlabel('Exoplanet')
    resq1.set_ylabel(r'$q_1$ c - u')
    
    # PLOT 2
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(4, 2, height_ratios=[2,2,1,1])
    
    # q0
    q0 = fig.add_subplot(gs[0])
    q0.scatter(q0_coupled, q0_uncoupled, marker='.', lw=0, alpha=0.5, edgecolor='none')
    q0.scatter(q0_coupled_top10, q0_uncoupled_top10, marker='.', lw=0, color='k')
    q0.set_xlabel(r'$q_0$ coupled')
    q0.set_ylabel(r'$q_0$ uncoupled')
    
    r = np.corrcoef(q0_coupled_top10, q0_uncoupled_top10)[0,1]
    from matplotlib.offsetbox import AnchoredText
    txt = AnchoredText(f'$R=${r:.2f}',
                       frameon=False,loc='lower right',
                       prop=dict(fontweight="bold"))
    q0.add_artist(txt)
    m, b = np.polyfit(q0_coupled_top10, q0_uncoupled_top10, 1)
    q0.plot(q0_coupled_top10, m*q0_coupled_top10+b, marker='', lw=0.5, color='k')
    # q1
    q1 = fig.add_subplot(gs[2])
    q1.scatter(q1_coupled, q1_uncoupled, marker='.', lw=0, alpha=0.5, edgecolor='none')
    q1.scatter(q1_coupled_top10, q1_uncoupled_top10, marker='.', lw=0, color='k')
    q1.set_xlabel(r'$q_1$ coupled')
    q1.set_ylabel(r'$q_1$ uncoupled')
    
    r = np.corrcoef(q1_coupled_top10, q1_uncoupled_top10)[0,1]
    txt = AnchoredText(f'$R=${r:.2f}',
                       frameon=False,loc='lower right',
                       prop=dict(fontweight="bold"))
    q1.add_artist(txt)
    m, b = np.polyfit(q1_coupled_top10, q1_uncoupled_top10, 1)
    q1.plot(q1_coupled_top10, m*q1_coupled_top10+b, marker='', lw=0.5, color='k')
    



