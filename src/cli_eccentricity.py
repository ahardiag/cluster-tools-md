#!/usr/bin/env python
# coding: utf-8

import os, sys, time, glob, re ,math
import argparse
from argparse import RawTextHelpFormatter
import warnings

import MDAnalysis as mda
import numpy as np
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
import pandas as pd
import scipy
from numba import jit
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)
from tqdm import tqdm
from time import sleep

def parse_args(required=False):
    parser = argparse.ArgumentParser(description="Compute Spherical Radial Density around cluster COM.",
                                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-s", "--sim",
                        required=required,
                        help="Name of the Data trajectory.")
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
    parser.add_argument('--method',
                        choices=['double_loop', 'best_centering'],
                        default='double_loop',
                        help='Choose the method for computing the inertia tensor.')
    parser.add_argument('--compound',
                        choices=['residues', 'atoms'],
                        default='residues',
                        help='Choose the level at which is computed the inertia tensor.\n'
                             '`residues` is faster as it only consider the center of mass of the molecules.') 
    parser.add_argument("--rename",
                        default=False,
                        help="Set to True to consider sim name with no suffix `Data_`.")
    parser.add_argument("--in_memory",
                        action="store_true",
                        help="Charge the whole trajectory in RAM memory.\n"
                        "Faster than reading from file.")
    parser.add_argument("-z", "--size",
                        type=int,nargs='+',
                        help="Choose specific size(s) of clusters.")
    parser.add_argument("--range",
                        action='store_true',
                        help="Range of values provided by three values `start end step`"
                        "by -z/--size option.")
    parser.add_argument("-v","--verbose",
                        action='count',
                        default=0,
                        help="Turn on for more details during the \n"
                        "calculation.")
    return parser.parse_args()

def ar2str(ar):
    string = ''
    for el in ar:
        string+=' '+str(el)
    return string 

def create_dir(args):
    os.makedirs(args.outdirname,exist_ok=True)
    if args.outfigdirname is None :
        args.outfigdirname = 'figures/{}'.format(os.path.basename(os.path.normpath(args.csvdir)))
    os.makedirs(args.outfigdirname,exist_ok=True)

@jit(nopython=True)
def delta(r1,r2,box=None,method='euclidian'):
    """
    Algebric distance between two points. Periodic boundary conditions can be settled
    by setting `method` keyword to *pbc*.
    Only orthogonal box is implemented.
    """
    if method == 'euclidian':
        return r2 - r1
    elif method == 'pbc':
        assert(box.shape == (3,)); "Box dimensions must be given in a 3 numpy array."
        dr   = r2 - r1
        sign = np.sign(dr)
        dist = sign*dr
        cond = (dist > box / 2)
        return sign*(dist - box)*cond + dr*~cond 
    else :
        raise ValueError("The method for calculating distances can only be *euclidian* or *pbc*.")

@jit(nopython=True)
def inertia_double_loop(coordinates, masses, total_mass=None,box=None,method=None):
    """
    Compute the inertia tensor using the distances between atoms and a double loop implementation.
    
    The advantage of this method is that it does not refer to the center of mass of the aggregate 
    which can be difficult to find when the aggregate is divided by the periodic box.
    """
    # coordinates change for each frame
    n_atoms = coordinates.shape[0]
    S       = np.zeros(6)
    prod = np.zeros(6)
    # double loop with non repeating pairs
    for i in range(n_atoms):
        for j in range(i+1,n_atoms) :
            # algebraic distances
            drij = delta(coordinates[i],coordinates[j],box=box,method=method) 
            # diagonal elements
            prod[0] = drij[0]*drij[0]
            prod[3] = drij[1]*drij[1] 
            prod[5] = drij[2]*drij[2] 
            sumprod = prod[0] + prod[3] + prod[5]
            prod[0] = sumprod - prod[0] # xx
            prod[3] = sumprod - prod[3] # xy
            prod[5] = sumprod - prod[5] # zz
            # products elements
            prod[1] = - drij[0]*drij[1] # xy
            prod[2] = - drij[0]*drij[2] # xz
            prod[4] = - drij[1]*drij[2] # yz
            S += prod / (total_mass) * masses[i] * masses[j] 
    return S

def radgyr_inertia_dist(atomgroup,method=None,debug=False,compound='atoms'):
    """
    Compute the total radius of gyration and the principal moments using the inertia tensor obtained
    with distances.
    
    Numpy arrays are required as inputs.
    
    Parameters
    ----------
    method : {'euclidian', 'pbc'}
        If ``'euclidian'``, the distance between positions do not take into account the periodic box.
        Otherwise, take into account the periodic box.
    compound : {'atoms', 'segments', 'residues', 'molecules', 'fragments'},\
           optional
    If ``'atoms'``, the positions of all atoms in the group will
    be used as a single position vector. Otherwise, the centers of
    mass of each :class:`Segment`, :class:`Residue`, molecule, or
    fragment will be used as an array of position vectors, i.e. a 2d
    array.    
    
    See also : inertia_double_loop function
    """
    # Reduce the dimensionality
    if compound == "atoms":
        coordinates, masses, total_mass  = atomgroup.positions, atomgroup.masses, atomgroup.masses.sum()
    else :
        group = getattr(atomgroup,compound)
        coordinates  = atomgroup.center_of_mass(compound=compound)
        masses, total_mass = group.masses, group.masses.sum()
    box = atomgroup.dimensions[:3] # TO DO : triclinic box
    
    I  = squareform_diagfill(inertia_double_loop(coordinates=coordinates,masses=masses,total_mass=total_mass,
                                                 box=box,method=method))
    values, evecs = np.linalg.eigh(I)
    indices = np.argsort(values)
    U = evecs[:, indices]
    Lambda = U.T.dot(I.dot(U))
    
    if debug :
        print("Inertia Tensor in the diagonal basis : \n{}".format(Lambda))
        print("Differences : \n{}".format(Lambda - np.diag(np.diagonal(Lambda))))
    assert np.allclose(Lambda - np.diag(np.diagonal(Lambda)),0,atol=1e-07), "Problem of diagonalization of inertia tensor!"
    
    l1_square,l2_square,l3_square = (Lambda[i,i] for i in range(3))
    if (l1_square < 0)|(l2_square < 0)|(l3_square < 0):
        return np.array([np.nan,np.nan,np.nan,np.nan])
    else :
        rg = np.sqrt((l1_square +  l2_square + l3_square) / 2 / total_mass)
        return np.array([rg,np.sqrt(l1_square)/total_mass,np.sqrt(l2_square)/total_mass,np.sqrt(l3_square)/total_mass])

def squareform_diagfill(arr1D):
    #https://stackoverflow.com/questions/38747311/how-to-create-a-symmetric-matrix-from-a-numpy-1d-array-the-most-efficient-way
    n = int(np.sqrt(arr1D.size*2))
    if (n*(n+1))//2!=arr1D.size:
        print("Size of 1D array not suitable for creating a symmetric 2D array!")
        return None
    else:
        R,C = np.triu_indices(n)
        out = np.zeros((n,n),dtype=arr1D.dtype)
        out[R,C] = arr1D
        out[C,R] = arr1D
    return out

def eccentricity(rg_array):
    """
    Compute eccentricity which measures the anisotropy of the mass distribution of an aggregate giving 
    its principal moments.
    
    Tends to 0 for an sperical object.
    Tends to 1 for a 1D object.
    """
    rg,rx,ry,rz = rg_array
    r2 = rx*rx + ry*ry + rz*rz
    if r2 == 0:
        return np.nan
    else:
        return 1 - min(rx*rx,ry*ry,rz*rz)/(1/3 * r2)

def plot_ecc(df,filename,norm_x=None,
             logx=True,
             marker='^',ms = 6,alpha_err=0.2,
             left=0.2,bottom=0.3,right=0.8,top=None):

    fig,ax = plt.subplots()
        
    x   = df.to_numpy()[:,0]
    y   = df.to_numpy()[:,1]
    std = df.to_numpy()[:,2]
    
    if norm_x is not None : x = x/norm_x
    markers, caps, bars = ax.errorbar(x, y, yerr=std,ms=ms,capsize=5,
                capthick=2,fillstyle='none',markeredgewidth=2)

    # loop through bars and caps and set the alpha value
    [bar.set_alpha(alpha_err) for bar in bars]
    [cap.set_alpha(alpha_err) for cap in caps]

    plt.xlabel('n/N') if norm_x else plt.xlabel('aggregation size')
    plt.ylabel("eccentricity")
    if norm_x is not None: plt.xlim(0,1)
    if logx: plt.xscale('log')
    plt.subplots_adjust(left=left,bottom=bottom,right=right,top=top)
    plt.show()
    fig.savefig(filename,facecolor="w",dpi=200)
    return fig

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

    def localize_cluster_all_sizes(self,time_begin=0,time_end=500):
        '''
        Get a list of a 2-list (time and resids) of a cluster aggregate for a selected size.

        One can restrict the time interval where to search by using `time_begin` and `time_end` (in nanoseconds).
        '''
        df_sub = self.df_resids.loc[self.sim]
        
        # get the list of possible sizes
        sizes = df_sub.index
        
        # Check that time intervals matchs with the data in Dataframe with clustering results
        #assert time_begin > float(self.df_resids.columns[0]) ,'Parameters `time_begin` must be larger than the smaller time in columns of the input CSV file.'
        #print(float(self.df_resids.columns[-1]))
        #assert time_end < float(self.df_resids.columns[-1]) ,'Parameters `time_end` must be smaller than the smaller time in columns of the input CSV file.'
        
        loc_cluster = {size:[] for size in sizes}
        
        for time in self.df_resids.columns:
#            if (float(time)<time_begin)|(float(time)>time_end):
#                continue
            for label in np.unique(df_sub[time]):
                test = df_sub[time]==label
                loc_cluster[sum(test)].append([time,list(df_sub[time].index[test])])
        
        # erase all empty keys
        for size in sizes:
            if len(loc_cluster[size])==0:
                del(loc_cluster[size])
        return loc_cluster

    def get_eccentricity_double_loop(self,ag,nmol=None,
                                    selection_solute="resname H01",
                                    time_begin=0,time_end=500,
                                    frame_per_ns=1,verbose=True,
                                    compound='atoms',debug=False):
        """
        Compute the eccentricity of an aggregate for a single trajectory.
        
        Parameters
        ----------
        nmol : int or array or tuple or list or None, optional
            Number of molecules of solute.
        selection_solute : str
            Selection in MDAnalysis format for the solute atoms.
        time_begin : float
            Initial time to perform analysis.
        time_end : float
            Final time to perform analysis.    
        frame_per_ns : int, optional
            Number of frames per nanosecond of simulation (sampling frequency) 
        verbose : bool, optional
            If True, print the sizes of the aggregates during analysis
        compound : {'atoms', 'segments', 'residues', 'molecules', 'fragments'},\
            optional
        If ``'atoms'``, the positions of all atoms in the group will
        be used as a single position vector. Otherwise, the centers of
        mass of each :class:`Segment`, :class:`Residue`, molecule, or
        fragment will be used as an array of position vectors, i.e. a 2d
        array.

        returns
        ----------
        ecc : numpy array with shape (n_sizes,3) with n_sizes the total number of sizes detected in all frames.
            The second dimension is organized as follows :
            0: size of aggregates
            1: mean value of eccentricity 
            2: standard deviation

        """
        if type(nmol).__name__ == 'NoneType':
            # Guess number of solute molecules
            nmol = np.unique(ag.select_atoms(selection_solute).resids).shape[0]
            sizes = range(1,nmol+1)
        elif type(nmol).__name__ == 'int':
            sizes = range(1,nmol+1)
        elif type(nmol).__name__ == 'str':
            raise TypeError; "Input `nmol` must be None, integer, tuple, list or array-like. "
        else:
            sizes = tuple(nmol)
        
        lsizes = len(list(sizes))
        if debug:
            print(nmol,lsizes)
            return
    
        ecc = []
        print('Reading the CSV input and find locators...',end='')
        locs = self.localize_cluster_all_sizes(time_begin=time_begin,time_end=time_end)
        print('Done.')
        
        # Compute only for the size specified by `nmol` parameter. 
        locs = {size:locs[size] for size in sizes if size in locs}
        
        #pb = progress_bar()
        pbar = tqdm(total=len(locs))
        for i,(size,loc) in enumerate(locs.items()):
            #loc = localize_cluster(sim,df_resids,size,time_begin=time_begin,time_end=time_end)
            if len(loc)==0: 
                continue
            #ecc[size] = []
            ecc_stats=[]
            for time,resids in loc:
                frame = int(float(time)*frame_per_ns) # to be improved to be not dependent on the way df_resids have been written
                ag_solute = ag.select_atoms("{} and resid {}".format(selection_solute,ar2str(resids)))
                ag_solute.universe.trajectory[frame]
                box = ag_solute.universe.dimensions[:3]
                
                ecc_stats.append(eccentricity(radgyr_inertia_dist(ag_solute,
                                                                method = 'pbc', compound=compound)))
            
            # Compute mean and standard deviation for a given aggregation size
            ecc.append([size,np.mean(ecc_stats),np.std(ecc_stats)])
            pbar.update(1)
        pbar.close()
        
        return np.array(ecc)


    def fill_df_ecc(self,args,
                  nmol=None,verbose=True,
                  debug=False,
                  time_begin=0,time_end=None,
                  method='double_loop',compound='atoms'):
        """
        Create a Dataframe with eccentricity values.
        
        Inputs 
        ----------
        args : namespace
            Arguments form argparse
        
        Parameters
        ----------
        nmol : int or array or tuple or list or None, optional
            Number of molecules of solute.
        method : {'double_loop', 'best_centering'}
            If ``'double loop'``, the calculation of the inertia tensor is done using double loop 
            and Numba accelerator. It allows to consider properly the periodic boundary conditions.
            If ``'best_centering'``, the calculation of the inertia tensor is done using MDAnalysis 
            tools. The algorithm use a preprocessing of the atom coordinates to find the best centering
            conditions which minimize the total radius of gyration. It does not take rigourously pbc but
            rather try to center the objects in the center of the box to avoid pbc crossings.
        compound : {'atoms', 'segments', 'residues', 'molecules', 'fragments'},\
            optional
            To be used with ``'double loop'``  method.
            If ``'atoms'``, the positions of all atoms in the group will
            be used as a single position vector. Otherwise, the centers of
            mass of each :class:`Segment`, :class:`Residue`, molecule, or
            fragment will be used as an array of position vectors, i.e. a 2d
            array.
        debug : bool
            Compute the eccentricity but do not fill the dictionary provided in input.
        time_begin : float
            Initial time to perform analysis.
        time_end : float
            Final time to perform analysis.    


        Returns
        ----------
        ecc_dict : dictionary with keys the keyword for each simualtion and with values as a Dataframe with columns:
            - <size>: size of aggregates
            - <mean>: mean value of eccentricity 
            - <std> : standard deviation

        """
            
        if method == "best_centering":
            in_memory=True
        elif method == "double_loop":
            in_memory=False
        else:
            raise NameError('Keyword `method` only take arguments in {double_loop,best_centering}')
            
        u = self.load_data_one(args)
        #u = load_data_one(sim,fileref=FILEREF,filetrj=FILETRJ,verbose=True,in_memory=in_memory)
        print("n_atoms = {} ; n_frames = {}".format(u.trajectory.n_atoms,u.trajectory.n_frames))

        #print(u.trajectory.n_atoms)
        frame_per_ns = 1000/u.trajectory.dt
        if time_end is None: 
            time_end = u.trajectory.n_frames/frame_per_ns
        #if nmol is None : nmol = convert_par(sim)['N'] # convert_par(sim)["N"]
        ag = u.select_atoms("resname H01")
        columns = ["size","mean","std"]
        if method == "best_centering":
            ecc = self.get_eccentricity_over_sizes(ag,nmol=nmol,
                                                    frame_per_ns=frame_per_ns,
                                                    verbose=True)
        elif method == "double_loop":
            ecc = self.get_eccentricity_double_loop(ag,nmol=nmol,
                                                    time_begin=time_begin,time_end=time_end,
                                                    frame_per_ns=frame_per_ns,
                                                        verbose=True,compound=compound)

        return pd.DataFrame(ecc,columns=columns)

def main():
    args = parse_args(required=True)

    # Create output directories
    create_dir(args)

    simulation = Simulation(args.sim)
    
    # Get the dataframe with the infos on aggregates
    simulation.read_csv(args)
    
    # Load trajectory
    #u = simulation.load_data_one(args)

    # Consider a subset of aggregate size
    if args.range: 
        assert (len(args.size)==2)|(len(args.size)==3) , '--size option accept only two or three values.'
        if (len(args.size)==2):
            sizes = np.arange(args.size[0],args.size[1])
        else :
            sizes = np.arange(args.size[0],args.size[1],args.size[2])
    else :
        sizes = args.size
    
    # Build a dictionary indexed by the simulation key
    # TODO : implement sequential analysis of several trajectories
    df = simulation.fill_df_ecc(args,nmol=sizes,verbose=True,
                                method=args.method,compound=args.compound,
                                time_begin=None,time_end=None,debug=False) 
    
    plot_ecc(df,'{}/eccentricity_{}.png'.format(args.outfigdirname,args.sim))

if __name__ == "__main__":
    main()
