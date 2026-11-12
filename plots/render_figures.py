"""Render MIMIC figures from saved pooled estimates; no models are fitted."""
from pathlib import Path
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(os.environ.get('MIMIC_CODE_ROOT', Path(__file__).resolve().parents[1]))
OUT = Path(os.environ.get('MIMIC_OUTPUT_DIR', ROOT/'outputs'))
DATA = Path(os.environ.get('MIMIC_DATA_DIR', ROOT/'data'))
FIG = OUT/'figures'
COLORS = ['#0072B2', '#D55E00']
STYLES = ['-', (0, (4, 2.5))]
LABELS = ['No reduction identified', 'Reduced']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
    'pdf.fonttype':42,'svg.fonttype':'none'})

def read(path):
    return pd.read_csv(OUT/path)

def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    for extension in ['pdf', 'svg', 'png']:
        fig.savefig(FIG/f'{name}.{extension}',dpi=600,bbox_inches='tight')
    plt.close(fig)

def legend(ax):
    handles=[Line2D([0],[0],color=COLORS[g],ls=STYLES[g],label=LABELS[g]) for g in (0,1)]
    ax.legend(handles=handles,title='RV systolic function',frameon=False,fontsize=8)

def curve(ax, frame, x='pco2', value='estimate', scale=100):
    for g in (0,1):
        f=frame[frame.rv.eq(g)].sort_values(x)
        ax.fill_between(f[x],scale*f.lower,scale*f.upper,color=COLORS[g],alpha=.14,lw=0)
        ax.plot(f[x],scale*f[value],color=COLORS[g],ls=STYLES[g],lw=1.6)
    ax.grid(axis='y',color='#e5e5e5',lw=.5)

def rcs():
    d=read('contrasts/tables/mimic_risks.csv')
    r=read('contrasts/tables/mimic_relative40.csv')
    tests=read('contrasts/tables/mimic_tests.csv').set_index('test').p_value
    fig,axes=plt.subplots(2,1,figsize=(6.6,6.5),sharex=True,layout='constrained')
    curve(axes[0],d);curve(axes[1],r)
    axes[0].set(title='A   Standardized 28-day mortality',ylabel='Mortality (%)',ylim=(0,100))
    axes[1].set(title='B   Risk difference relative to 40 mmHg',ylabel='Risk difference (percentage points)',xlabel='PaCO₂ (mmHg)',xlim=(25,65))
    axes[1].axhline(0,color='#777777',ls='--',lw=.7)
    legend(axes[0])
    axes[0].text(.02,.02,f"Global interaction P = {tests['interaction']:.5f}\n"
        f"Nonlinearity P: no reduction = {tests['rv0_nonlinear']:.3f}; reduced = {tests['rv1_nonlinear']:.3f}",
        transform=axes[0].transAxes,fontsize=8,va='bottom')
    nadir_file=OUT/'contrasts/tables/RV0_minimum_risk_points.csv'
    if nadir_file.exists():
        nadir=float(pd.read_csv(nadir_file).paco2_minimum.iloc[0])
        for ax in axes:
            ax.axvline(nadir,color=COLORS[0],ls='--',lw=.8)
            ax.text(nadir,.02,f'{nadir:.1f}',transform=ax.get_xaxis_transform(),color=COLORS[0],ha='center',fontsize=8)
    save(fig,'rcs_standardized_mortality')

def interaction():
    d=read('contrasts/tables/mimic_sliding.csv')
    fig,axes=plt.subplots(2,1,figsize=(6.6,6.4),sharex=True,layout='constrained')
    for g in (0,1):
        f=d[d.metric.eq(f'change_rv{g}')].sort_values('x_start')
        axes[0].plot(f.x_start,100*f.estimate,color=COLORS[g],ls=STYLES[g])
        axes[0].fill_between(f.x_start,100*f.lower,100*f.upper,color=COLORS[g],alpha=.14,lw=0)
    f=d[d.metric.eq('interaction')].sort_values('x_start')
    axes[1].plot(f.x_start,100*f.estimate,color='#172A36')
    if 'simultaneous_lower' in f:
        axes[1].fill_between(f.x_start,100*f.simultaneous_lower,100*f.simultaneous_upper,color='#172A36',alpha=.08,label='Approximate simultaneous 95% band')
    axes[1].fill_between(f.x_start,100*f.lower,100*f.upper,color='#172A36',alpha=.20,label='Pointwise 95% interval')
    for ax in axes:
        ax.axhline(0,color='#777777',ls='--',lw=.7)
        ax.grid(axis='y',color='#e5e5e5',lw=.5)
    axes[0].set(title='A   Risk change from x to x + 10 mmHg',ylabel='Risk change (percentage points)')
    axes[1].set(title='B   Difference between RV groups',ylabel='Reduced minus no reduction\n(percentage points)',xlabel='Starting PaCO₂ (mmHg)',xlim=(25,55))
    legend(axes[0]);axes[1].legend(frameon=False,fontsize=8)
    save(fig,'continuous_additive_interaction')

def sensitivity():
    fig,axes=plt.subplots(1,2,figsize=(8.0,3.8),sharex=True,sharey=True,layout='constrained')
    for ax,folder,tag,label in zip(axes,['contrasts','sensitivity_6h'],['mimic','mimic6h'],['Main: ±12 h','Sensitivity: ±6 h']):
        d=read(f'{folder}/tables/{tag}_risks.csv');p=read(f'{folder}/tables/{tag}_tests.csv').set_index('test').loc['interaction','p_value']
        curve(ax,d);ax.set(title=label,xlim=(25,65),ylim=(0,100),xlabel='PaCO₂ (mmHg)')
        ax.text(.03,.02,f'Global interaction P = {p:.5f}',transform=ax.transAxes,fontsize=8)
    axes[0].set_ylabel('Standardized 28-day mortality (%)');legend(axes[0])
    save(fig,'pairing_window_sensitivity')

def support():
    d=read('support/tables/forest.csv')
    fig,axes=plt.subplots(1,2,figsize=(10.5,6.0),sharex=True,layout='constrained')
    for ax,days in zip(axes,[7,28]):
        labels=[]
        for row,(stem,label) in enumerate([('vpfd','Vasopressors'),('vfd','Ventilation'),('crrtfd','CRRT')]):
            f=d[d.outcome.eq(f'{stem}{days}')]
            p=float(f.interaction_p.iloc[0])
            for k,category in enumerate(['low','high']):
                for g in (0,1):
                    r=f[f.category.eq(category)&f.rv.eq(g)].iloc[0]
                    y=2*row+k+(-.12 if g==0 else .12)
                    ax.errorbar(r.estimate,y,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],
                        fmt='o' if g else 's',color=COLORS[g],capsize=2,ms=4)
            labels.extend([f'{label}: <35\nInteraction P = {p:.3f}',f'{label}: >45'])
        ax.set_yticks(range(6),labels,fontsize=8)
        ax.set_ylim(5.6,-.6);ax.set_title(f'{days}-day support-free days')
        ax.axvline(0,color='#777777',ls='--',lw=.7)
        ax.set_xlabel('Adjusted difference versus 35–45 mmHg (days)')
    handles=[Line2D([0],[0],marker='s' if g==0 else 'o',ls='',color=COLORS[g],label=LABELS[g]) for g in (0,1)]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1.07),ncol=2,frameon=False,title='RV systolic function')
    save(fig,'support_free_days')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('figure',choices=['rcs','interaction','sensitivity','support','all'])
    selection=parser.parse_args().figure
    functions={'rcs':rcs,'interaction':interaction,'sensitivity':sensitivity,'support':support}
    for name in functions if selection=='all' else [selection]:
        functions[name]()
