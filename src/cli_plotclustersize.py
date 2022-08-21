#!/usr/bin/env python 

import MDAnalysis as mda
import numpy as np
import os, sys, time, glob, re ,math

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl
#import seaborn as sns
import pandas as pd

import argparse
from argparse import RawTextHelpFormatter

def parse_args(required=True):
    parser = argparse.ArgumentParser(description="Test plots from CSV files obtained with clustsize program.",
                                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--inputdir",
                        default="outputs",
                        help="Directory with CSV files.")
    parser.add_argument("-o", "--outputdir",
                        default="figures",
                        help="Directory where to create figures.")
    parser.add_argument("--tmin",
                        default=0,
                        type=float,
                        help="Start time for timeseries (in ns)")
    parser.add_argument("--tmax",
                        required=required,
                        type=float,
                        help="Start time for timeseries (in ns)")
    parser.add_argument("-v","--verbose",
                        action='store_true',
                        help="Turn on for more details during the calculation.")
    parser.add_argument("--rename",
                        action='store_true',
                        help="Option to deal with simulation names saved with sufffix <Data_>")
    return parser.parse_args()

def create_dir(args):
    try:
        os.makedirs(args.outputdir)
        if args.verbose: print("Directory " , args.outputdir ,  " Created ")
    except FileExistsError:
        if args.verbose: print("Directory " , args.outputdir ,  " already exists")  
 
def get_formatted_name(sim,rename):
    '''
    Get rid of the suffix if data is located in a directory starting with <Data_> string. 
    '''
    if rename:
        if sim[0:5]=="Data_":
            sim=sim[5:]
    return sim

def convert_par(sim,rename):
    '''
    Convert the name of the data directory into the set of pertinent parameters for the analysis.
    '''
    sim = get_formatted_name(sim,rename=rename)
    N   = int(sim.split("_")[1])
    m   = sim.split("_")[3].split("m")[1]
    box = float(sim.split("_")[4].split("box")[1])
    if (m=='inf') :
        f = 1
    elif (m=='0'):
        f = 0
    else:
        f   = 1/(1+1/float(m))
    C   = float(N)/(0.0006022*(box)**3)
    keys=["N","m","box","C","f"]
    values=[N,m,box,C,f]
    return dict(zip(keys,values))

def get_dict_colors(sim_name,rename):
    #colors={get_formatted_name(sim,rename):sns.color_palette('tab20',n_colors=len(sim_name))[ic] for ic,sim in enumerate(sim_name)}
    colors={get_formatted_name(sim,rename):plt.get_cmap('tab20')(ic) for ic,sim in enumerate(sim_name)}
    return colors

def extract_sim_names(args):
    df_prop=pd.read_csv('{}/cluster_properties.csv'.format(args.inputdir))
    df_prop.set_index(['simulation','property'],inplace=True)
    return df_prop.index.levels[0]

def read_csv_prop(args):
    df_prop=pd.read_csv('{}/cluster_properties.csv'.format(args.inputdir))
    df_prop.set_index(['simulation','property'],inplace=True)
    df_prop.columns=df_prop.columns.astype('float')
    df_prop.T.index.name='time (ns)'
    return df_prop

def read_csv_histo(args):
    df_histo=pd.read_csv('{}/cluster_histo.csv'.format(args.inputdir))
    df_histo.set_index(['simulation','property'],inplace=True)
    df_histo.shape,df_histo.columns.dtype
    df_histo=df_histo.T
    return df_histo

def plot_maxclust(df_prop,args,sim_name):
    prop='maxclust'
    plt.figure(figsize=(4,3))
    plt.subplots_adjust(bottom=0.3,left=0.2,right=0.8)
    list_sim=[sim for sim in sim_name if (sim.find('100')>0)]
    list_sim=sorted(list_sim, key=lambda sim: convert_par(sim,args.rename)['C'])

    label_sim=["%2.0f"%convert_par(sim,args.rename)["C"] for sim in list_sim]
    colors = get_dict_colors(sim_name,args.rename)

    ymax=0
    for sim in list_sim:
        x=df_prop.T.loc[:,([sim],prop)].index.to_numpy()
        y=df_prop.T.loc[:,([sim],prop)].to_numpy()
        plt.plot(x,y,color=colors[get_formatted_name(sim,args.rename)])
        ymax=max(ymax,y.max())
    plt.legend(label_sim,ncol=3,fontsize=6,title='Concentration (mmol/L)')
    plt.xlabel("time (ns)")
    plt.ylabel("Normalized Size \n of the biggest cluster",fontsize=8)
    plt.xlim(args.tmin,args.tmax)
    plt.ylim(-0.2,ymax*1.01)
    filename = '{}/maxclust-f{:0.2f}.png'.format(args.outputdir,convert_par(sim,args.rename)["f"])
    plt.savefig(filename,dpi=200)
    print("File %s created !"%filename)

def plot_histo(df_histo,args,sim_name):
    plt.figure(figsize=(4,3))
    list_sim=[sim for sim in sim_name if (sim.find('100')>0)]
    list_sim=sorted(list_sim, key=lambda sim: convert_par(sim,args.rename)['C'])

    label_sim=["%2.0f"%convert_par(sim,args.rename)["C"] for sim in list_sim]
    colors = get_dict_colors(sim_name,args.rename)
    for sim in list_sim:
        x=df_histo.loc[:,sim].to_numpy()[:,0]
        y=df_histo.loc[:,sim].to_numpy()[:,1]
        plt.plot(x,y,color=colors[get_formatted_name(sim,args.rename)],marker="^",lw=0,markersize=5)
    plt.yscale('log')
    plt.ylabel('P(n)')
    plt.xlabel('n')
    plt.legend(label_sim,ncol=3,fontsize=7,title='Concentration (mmol/L)')
    plt.subplots_adjust(left=0.2,bottom=0.3,right=0.8)
    #plt.show()
    filename = '{}/histo-clust-f{:0.2f}.png'.format(args.outputdir,convert_par(sim_name[0],args.rename)["f"])
    plt.savefig(filename,dpi=200)
    print("File %s created !"%filename)

def __init__(self):
    return None

def main():
    args = parse_args(required=True)
    create_dir(args)
    sim_name = extract_sim_names(args)
    #print(sim_name)
    #sys.exit()
    df_prop = read_csv_prop(args)
    df_histo = read_csv_histo(args)

    # Plot the time serie of the biggest cluster 
    plot_maxclust(df_prop,args,sim_name)

    # Plot the histogram distrbution of cluster size in the chosen time window 
    plot_histo(df_histo,args,sim_name)

if __name__ == "__main__":
    main()