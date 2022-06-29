#!/usr/bin/env python
# coding: utf-8
import os, sys, time, glob, re ,math
import gc 

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style

import numpy as np
import seaborn as sns
import pandas as pd
import scipy
import MDAnalysis as mda

import argparse
from argparse import RawTextHelpFormatter

#https://github.com/tqdm/tqdm/
from tqdm import tqdm
from time import sleep

from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)

#assert int(mda.__version__[0])>1 , 'Script need version of MDAnalysis > 2.0.0'

def parse_args(required=False):
    parser = argparse.ArgumentParser(description="Compute Spherical Radial Density around cluster COM.",
                                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-s", "--sim",
                        required=required,
                        help="Name of the Data trajectory.")
    parser.add_argument("-z", "--size",
                        required=required,
                        type=int,nargs='+',
                        help="Size of the cluster.")
    parser.add_argument("--range",
                        action='store_true',
                        help="Range of values provided by three values `start end step`"
                        "by -z/--size option.")
    parser.add_argument("-i", "--inputdir",
                        required=required,
                        help="directory with all trajectories to analyse.")
    parser.add_argument("-c", "--csvdir",
                        required=required,
                        help="directory with the dataframe `cluster_resids.csv`.")
    parser.add_argument("--nbins",
                        type=int,
                        default=50,
                        help="Number of bins of the density distribution.")
    parser.add_argument("-m", "--resname",
                        default="H01",
                        help="Resname or molecule name. e.g: HC8.")
    parser.add_argument("-o", "--outdirname",
                        default="./outputs/",
                        help="Output Sub-directory name.")
    parser.add_argument("--outfigdirname",
                        help="Output Sub-directory name for figures.")
    parser.add_argument("-r", "--fileref",
                        default="init.gro",
                        help="Commmon name of reference file in all \n"
                        "simulation directories.")
    parser.add_argument("-f", "--filetrj",
                        default="md_trj.xtc",
                        help="Commmon name of trajectory file in all \n"
                        "simulation directories.")
    parser.add_argument("--rename",
                        default=False,
                        help="Set to True to consider sim name with no suffix `Data_`.")
    parser.add_argument("--in_memory",
                        action="store_true",
                        help="Charge the whole trajectory in RAM memory.\n"
                        "Faster than reading from file.")
    parser.add_argument("-v","--verbose",
                        action='count',
                        default=0,
                        help="Turn on for more details during the \n"
                        "calculation.")
    return parser.parse_args()

def ar2str(ar,sep=' '):
    string = ''
    for el in ar:
        string+=sep+str(el)
    return string 

def create_dir(args):
    os.makedirs(args.outdirname,exist_ok=True)
    if args.outfigdirname is None :
        args.outfigdirname = 'Figures/{}'.format(os.path.basename(os.path.normpath(args.csvdir)))
    os.makedirs(args.outfigdirname,exist_ok=True)
    
def get_resids_from_labels(labels):
    '''
    Return a list of array with resids for each label cluster.
    Input must be a from a (2,nmol) array with row 1 : resids  row 2: labels
    Inputs must be sorted along the second row.
    
    '''
    assert(labels.shape[0]==2), "Input array must be of shape (2,nmol)."
    assert(np.all(np.diff(labels[1]) >= 0)) , "INput array must be sorted along the second row"

    delimiter=np.where( np.diff(labels[1,:])!= 0)[0]
    ar_out=np.split(labels[0,:],delimiter+1)
    return ar_out

def get_selection_from_resids(ar_index,ag):
    '''
    Return MDAnalysis selection name and atomgroups from the list of array with resids for each cluster.
    The universe
    '''
    ag_list=[]
    selections=[]
    for ar in ar_index:
        string="resname H01 and resid "
        for el in list(ar):
            string+=str(el)+" "
        selections.append(string)
        ag_list.append(ag.select_atoms(string))
    return ag_list,selections

def get_labels_from_df(sim,df_resids,index_time):
    '''
    Returns labels from dataframe with correct sorting with respect to label ID row and time in ns.
    '''
    df_tmp=df_resids.loc[([sim]),:].iloc[:,index_time].reset_index('resid')
    time=df_tmp.columns[1]
    df_tmp.sort_values(by=time,inplace=True)
    return df_tmp.to_numpy().T,time

def get_selection_from_df(sim,df_resids,frame,ag):
    '''
    Returns list resids and atomgroups from dataframe and time in ns.
    '''
    labels,time=get_labels_from_df(sim,df_resids,frame)
    ar_index=get_resids_from_labels(labels)
    return get_selection_from_resids(ar_index,ag),time

def center_on_subgroup(ag,selection,box):
    center = ag.select_atoms(selection).atoms.center_of_mass()
    ag.translate(box/2 - center)
    ag.wrap(compound="residues")

def distance(a, b, dimensions):
    """
    Get displacement in each coordinate and wrap w.r.t. lattice parameter.
    """
    dx = abs(a[0] - b[0])
    x = min(dx, abs(dimensions[0] - dx))
     
    dy = abs(a[1] - b[1])
    y = min(dy, abs(dimensions[1] - dy))
     
    dz = abs(a[2] - b[2])
    z = min(dz, abs(dimensions[2] - dz))
    return np.sqrt(x**2 + y**2 + z**2)

def plot_dist(counts,fontsize=12,figsize=(5,3),
              norm=True,savefig=False,
              r_min=None,
              filename='density.png',**kwargs):
    if norm:
        counts[1] /= counts[1].max()
        counts[2] /= counts[2][-1]
        counts[3] /= counts[3][-1]
    plt.figure(figsize=figsize)
    if r_min:
        cond = counts[0]>r_min
        counts = counts[:,cond]
    plt.plot(counts[0],counts[1],label='solute',**kwargs)
    plt.plot(counts[0],counts[2],label='water',**kwargs)
    plt.plot(counts[0],counts[3],label='ethanol',**kwargs)
    plt.xlabel(r"$r(\AA)$",fontsize=fontsize)
    plt.ylabel("Relative Density",fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.tight_layout()
    if savefig:
        plt.savefig(filename,dpi=300,facecolor='white', edgecolor='white')
    plt.show()

def compute_radial_density(u,loc,resolution=20,com='residues',r_cutoff=None,compound='H01',verbose=False):
    '''
    Computes the radial density distribution of solute,ethanol and water around the COM of an aggregate.
    
    Returns a 2-list as an 1-array for radius distances and a ndarray with the radial
    density of 3 types of atoms (solute,ethanol,water)
    '''
    if r_cutoff is None:
        r_cutoff = u.trajectory.ts.dimensions[0]/2 # cutoff for maximal distance
    dt = u.trajectory.dt
    dr = r_cutoff / resolution

    # volume of a sphere slice
    vol = np.zeros(resolution)
    for j in range(resolution):
        vol[j] = (j+1)*(j+1)*dr*dr*dr

    counts = np.zeros([len(loc),4,resolution])
    
    pbar = tqdm(total=len(loc))
    for i,(time,resids) in enumerate(loc):
        frame = int(float(time)*1000/dt)
        if verbose : print("moving to frame {} at time {}".format(frame,time))
        u.trajectory[frame]
        box = u.trajectory.ts.dimensions[:3]
        
        # center on a random resid
        ag_all = u.select_atoms("all")
        center_on_subgroup(ag_all,'resname {} and resid {}'.format(compound,resids[0]),box)

        center = u.select_atoms('resname {} and resid {}'.format(compound,ar2str(resids))).atoms.center_of_mass()
        
        ag_solute  = u.select_atoms('resname {} and resid {}'.format(compound,ar2str(resids)))
        ag_ethanol = u.select_atoms('resname ETH')
        ag_water   = u.select_atoms('resname SOL TIP3')
        
        if com == 'residues':
            coords_solute  = np.array([res.atoms.center_of_mass() for res in ag_solute.residues])
            coords_ethanol = np.array([res.atoms.center_of_mass() for res in ag_ethanol.residues])
            coords_water   = np.array([res.atoms.center_of_mass() for res in ag_water.residues])
        elif com == 'atoms':
            coords_solute  = ag_solute.atoms.positions
            coords_ethanol = ag_ethanol.atoms.positions
            coords_water   = ag_water.atoms.positions
        else:
            raise NameError('`compound` option can only take two values : `atoms` and `residues` ')
        
        # Loop on all atom positions
        for imol,coords in enumerate([coords_solute,coords_water,coords_ethanol]):
            for position in coords:
                dist = distance(position, center,box)
                index = int(dist / dr)
                if 0 < index < resolution:
                    counts[i,imol+1,index] += 1.0
        pbar.update(1)
    pbar.close()

    # compute total number of molecules
    for imol,coords in enumerate([coords_solute,coords_water,coords_ethanol]):
        counts[:,imol+1] /= vol*(counts[:,imol+1].sum()/(r_cutoff**3))
    
    # distance array
    counts[:,0,:] = np.arange(0,r_cutoff,dr)
    return counts

class Simulation():

    def __init__(self,sim):
        """
        Set up the initial analysis parameters.
        """
        self.sim = sim
    
    def load_data_one(self,args,verbose=True):
        if verbose : print("\nLoading trajectory from parent {}".format(self.sim))
        fileref=args.inputdir + self.sim + "/" + args.fileref
        filetrj=args.inputdir + self.sim + "/" + args.filetrj

        u = mda.Universe(fileref,filetrj,in_memory=args.in_memory,tpr_resid_from_one=False)
        assert u.atoms.center_of_mass().shape[0]>0, 'No atoms in trajectory !\n'
        assert u.trajectory.n_frames > 1 , "No frames loaded in trajectory, wrong file name !?"
        return u
    
    def read_csv(self,args):
        # Read csv
        self.df_resids=pd.read_csv(args.csvdir+'/cluster_resids.csv')
        #df_resids=pd.read_csv(PATH_TO_CSV+'/cluster_resids_order.csv')
        self.df_resids.set_index(['simulation','resid'],inplace=True)
        self.df_resids.columns.names=["time (ns)"]
        print("Reading Dataframe with size {}".format(self.df_resids.shape))
        
    def localize_cluster(self,size):
        '''
        Get a list of a 2-list (frame and resids) of a cluster aggregate for the selected size 
        '''
        df_sub = self.df_resids.loc[self.sim]
        loc_cluster=[]
        for col in self.df_resids.columns:
            for label in np.unique(df_sub[col]):
                test = df_sub[col].to_numpy()==label
                if (test).sum()==size:
                    loc_cluster.append([col,df_sub[col].index[test].to_numpy()])
        return loc_cluster

def main():
    args = parse_args(required=True)

    # Create output directories
    create_dir(args)

    simulation = Simulation(args.sim)
    
    # Get the dataframe with the infos on aggregates
    simulation.read_csv(args)
    
    # Load trajectory
    u = simulation.load_data_one(args)
        
    # Loop on sizes
    if args.range: 
        assert (len(args.size)==2)|(len(args.size)==3) , '--size option accept only two or three values.'
        if (len(args.size)==2):
            sizes = np.arange(args.size[0],args.size[1])
        else :
            sizes = np.arange(args.size[0],args.size[1],args.size[2])
    else :
        sizes = args.size

    for size in sizes:
        # Get cluster frames from sizes
        loc = simulation.localize_cluster(size)
        if len(loc) == 0 :
            continue
        print("Aggregate of size {} : {} frames".format(size,len(loc)))
        
        # Compute the density around aggregate
        counts = compute_radial_density(u,loc,
                                        resolution=args.nbins,r_cutoff=None,
                                        com='residues',
                                        verbose=args.verbose)
        
        # Compute mean over all frames
        counts_mean = np.zeros([counts.shape[1],counts.shape[2]])
        counts_mean[0,:] = counts[0,0,:]
        counts_mean[1:,:] = counts[:,1:,:].mean(axis=0)
        
        # Plot the figure                                                    
        filename = '{}/density_{}_size{}.png'.format(args.outfigdirname,
                                                        args.sim,
                                                        size)
        plot_dist(counts_mean,norm=True,savefig=True,filename=filename,r_min=4)
        print("Plotting results in {}.".format(filename))

if __name__ == "__main__":
    main()
