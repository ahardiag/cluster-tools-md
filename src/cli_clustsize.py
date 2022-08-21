#!/usr/bin/env python 
import os, sys, time, glob, re ,math
import argparse
from argparse import RawTextHelpFormatter
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import self_distance_array,distance_array
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)
from sklearn.cluster import DBSCAN
from sklearn.neighbors import radius_neighbors_graph
from scipy.spatial.distance import squareform
import numpy as np
import pandas as pd
import fnmatch

class Clustsize(AnalysisBase):  # subclass AnalysisBase

    def __init__(self, atomgroup, verbose=True,debug=False,nmol=100,natoms=43,time_scale=1000,min_samples=5,eps=3,
                 with_distance_matrix=False,with_double_labelling=True):
        """
        Set up the initial analysis parameters.
        """
        # must first run AnalysisBase.__init__ and pass the trajectory
        trajectory = atomgroup.universe.trajectory
        super(Clustsize, self).__init__(trajectory,
                                               verbose=verbose)
        # set atomgroup as a property for access in other methods
        self.atomgroup = atomgroup
        # to compute only one time 
        self.n_frames = trajectory.n_frames
        self.nmol=nmol
        self.natoms=natoms
        self.ag_split_by_atom=atomgroup.split("atom")
        self.time_scale = time_scale
        self.min_samples = min_samples
        self.eps=eps = eps
        self.with_distance_matrix = with_distance_matrix
        self.n = self.nmol * self.natoms
        self.matrix_to_allocate = np.zeros(int(self.nmol * self.natoms * (self.nmol * self.natoms - 1) / 2))
        self.with_double_labelling = with_double_labelling
        self.verbose = verbose
        self.debug = debug

    def _prepare(self):
        """
        Create array of zeroes as a placeholder for results.
        This is run before we begin looping over the trajectory.
        """
        # This must go here, instead of __init__, because
        # it depends on the number of frames specified in run().
        
        # clustering statitistics (3,nmol) row:  1 for labels of clusters, 
        # 2 for number of atoms in clusters ,3 for resids list
        # columns: maximum = nmol if each
        # molecule is in a single cluster of size 1 (i.e. NATOMS data points)
        self.counts = np.zeros((self.n_frames, 3,self.nmol),dtype=int)
        
        # labels for atoms
        self.labels  = np.zeros((self.n_frames,self.nmol*self.natoms),dtype=int)
        self.labels_pbc  = np.zeros((self.n_frames,self.nmol*self.natoms),dtype=int)
        self.labels_new  = np.zeros((self.n_frames,self.nmol*self.natoms),dtype=int)
        
        # frames
        self.times = np.zeros(self.stop-self.start)

    def _single_frame(self):
        """
        This function is called for every frame that we choose
        in run().
        """
        # call our earlier function
        ac, frame_labels, frame_labels_pbc, frame_labels_new = analysis_cluster(self.atomgroup,self.ag_split_by_atom,
                                                                                nmol=self.nmol,natoms=self.natoms,
                                                                                min_samples=self.min_samples,eps=self.eps,
                                                                                with_distance_matrix=self.with_distance_matrix,
                                                                                matrix_to_allocate=self.matrix_to_allocate,
                                                                                with_double_labelling=self.with_double_labelling,
                                                                                verbose = self.verbose, debug = self.debug )
        
        # save labels and number of molecules in each cluster
        self.counts[self._frame_index] = ac
        # save labels of atoms
        self.labels[self._frame_index] = frame_labels
        self.labels_pbc[self._frame_index] = frame_labels_pbc
        self.labels_new[self._frame_index] = frame_labels_new
        
        # save times
        self.times[self._frame_index] = self._trajectory.time+1
        # save the list of resids for each cluster
        
    def _conclude(self):
        """
        Finish up by calculating an average and transforming our
        results into a DataFrame.
        """
        # by now self.result is fully populated
        #convert time 
        self.times /= self.time_scale 
        # save number of clusters
        self.nclust = np.count_nonzero(self.counts[:,1,:],axis=1)/self.nmol
        # save max size of clusters
        self.maxclust = np.max(self.counts[:,1,:],axis=1)/self.nmol
        # save histogram of sizes of clusters
        histo_data,bin_edges = np.histogram(self.counts[:,1,:],bins=self.nmol,range=[0.5,self.nmol+0.5])
        bin_centers = np.arange(1,self.nmol+1,1)
        self.histo_clust=np.asarray([bin_centers,histo_data])
        
        # check the normalisation of the data
        assert (histo_data*bin_centers).sum() == self.nmol * self.n_frames , "Normalization of the size histogram is wrong !"     

def parse_args(required=False):
    parser = argparse.ArgumentParser(description="Find all aggregate clusters for a set of MD trajectories.\n"
                                            "For selection commands, see MDAnalysis doc :\n"
                                            "https://docs.mdanalysis.org/stable/documentation_pagses/selections.html",
                                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--inputdir",
                        required=required,
                        help="directory with all trajectories to analyse")
    parser.add_argument("-s","--selection",
                        type=str,
                        help="MDAnalysis selection of solute atoms to apply the\n"
                        "distance clustering.Replace the -m option.")
    parser.add_argument("-m", "--resname",
                        default="H01",
                        help="Resname or molecule name.\n"
                        "e.g: HC8. Obsolete : use -s option instead")
    parser.add_argument("-o", "--outdirname",
                        default="./outputs/",
                        help="Output Sub-directory name.")
    parser.add_argument("-r", "--reffile",
                        default="init.gro",
                        help="Commmon name of reference file in all \n"
                        "simulation directories.")
    parser.add_argument("-f", "--trjfile",
                        default="md_trj.xtc",
                        help="Commmon name of trajectory file in all \n"
                        "simulation directories.")
    parser.add_argument("-b", "--begin",
                        type=int,
                        default=0,
                        help="Start frame of analysis.")
    parser.add_argument("-e", "--end",
                        required=required,
                        type=int,
                        help="Stop frame of analysis.")
    parser.add_argument("-v","--verbose",
                        action='count',
                        default=0,
                        help="Turn on for more details during the \n"
                        "calculation.")
    parser.add_argument("--mins",
                        required=required,
                        type=int,
                        help="Minimum number of points (atoms) to define\n"
                        " a cluster.")
    parser.add_argument("--eps",
                        required=required,
                        type=float,
                        help="Maximum distance between two points (atoms)\n"
                        "in the same cluster.")
    parser.add_argument("--with-distance-matrix",
                        action='store_true',
                        help="Use total matrix distance to treat periodic\n"
                        " boundary conditions. Slow for large systems!")
    parser.add_argument("--with-double-labelling",
                        action='store_true',
                        help="Use double labelling method to perform\n"
                        "clustering analysis on two shifted set of positions.\n"
                        "Allow to search for clusters close to the box boundaries.")
    parser.add_argument("--ignore-warnings",
                        action='store_true',
                        help="Avoid stdout all warnings from using\n"
                        "clustering on sparse matrix.")
    parser.add_argument("--pattern-dirs",
                        type=str,
                        default=None,
                        help="String pattern to select a serie of data\n"
                        "directories with similar pattern in the input directory.")    
    return parser.parse_args()

def check_version():
    assert int(mda.__version__[0])>1 , 'Script need version of MDAnalysis > 2.0.0'

def check_dirs(args):
    '''
    Detect the format of the directories with trajectory data.
    '''
    ldirs = os.listdir(args.inputdir)
    assert len(ldirs) >0 , 'No elements are found in the input directory.'
    lnum=[]
    
    if args.pattern_dirs is None:
        for dir in ldirs:
            lnum.append(dir.split('_')[0])
        assert np.unique(lnum).shape[0] == 1 ,\
        'All data directories must have the same structure with the same starting key (e.g : Data_...). '\
        'Use --pattern-dirs to use string pattern common in all directories to analyze.'
    else :
        ldirs = fnmatch.filter(ldirs,args.pattern_dirs)
        assert len(ldirs) >0 , 'No directory matches the following wildcard : "{}" !'.format(args.pattern_dirs)
    
    if args.verbose:
        print("Reading directories:\n",ldirs,flush=True)
    for sim in ldirs:
        assert len(glob.glob("{}/{}/{}".format(args.inputdir,sim,args.reffile)))>0, "No reference  file {}/{}/{} found !".format(args.inputdir,sim,args.reffile)
        assert len(glob.glob("{}/{}/{}".format(args.inputdir,sim,args.trjfile)))>0, "No trajectory file {}/{}/{} found !".format(args.inputdir,sim,args.trjfile)
    return ldirs

def check_method(args):
    if args.with_distance_matrix:
        print("\nMethod 1 : Total Matrix Distance. Compatible with periodic box. "
        "More accurate but slow since it needs to allocate full distance matrix "
        "~ o(N²).\n",flush=True)
    elif args.with_double_labelling:
        print("\nMethod 2 : Double Labelling. Compatible with periodic box. "
        "Two set of positions have been shifted by the half-diagonal box. "
        " Clustering is applied in both sets and the minimal common label "
        "is determined then.\n",flush=True)
    else :
        print("\nMethod 3 : Simple method . NOT COMPATIBLE WITH PERIODIC "
        "BOUNDARY CONDITIONS. Fast method since no distance matrix is allocated .\n",flush=True)

def check_warnings(args):
    if args.ignore_warnings:
        print("WARNING : Warnings will not appear in this choice of parameters !\n",flush=True)
        warnings.filterwarnings("ignore")
        
def create_dir(args):
    try:
        os.makedirs(args.outdirname)
        if args.verbose: print("Directory " , args.outdirname ,  " Created ",flush=True)
    except FileExistsError:
        if args.verbose: print("Directory " , args.outdirname ,  " already exists",flush=True) 
 
def load_traj(sim,args):
    fileref="{}/{}/{}".format(args.inputdir,sim,args.reffile)
    filetrj="{}/{}/{}".format(args.inputdir,sim,args.trjfile)
    if args.verbose: print("\n\nLoading trajectory files [{},{}] with MDanalysis :".format(args.reffile,args.trjfile),flush=True)
    u = mda.Universe(fileref,filetrj,in_memory=False)
    if args.verbose : 
        print("{:35s} {:8s} {:8s}".format("sim","n_frames","n_atoms"),flush=True)
        print("{:35s} {:8d} {:8d}".format(sim,u.trajectory.n_frames,u.trajectory.n_atoms),flush=True)
    assert u.atoms.center_of_mass().shape[0]>0, 'No atoms in trajectory !\n'
    return u

def shift_pbc(X,box_size):
    '''
    Substract the half box value of each dimension in the position coordinates.
    Only works if the data is centered on the center of the box.
    '''
    try:
        # Convert it into a float
        box_size = float(box_size)
        # Convert into an array
        box_size = np.array([box_size]*2)
    except TypeError:
        box_size= np.array(box_size)
        if box_size.shape[0] != X.shape[1]:
            print(box_size.shape[0],X.shape[1],flush=True)
            raise ValueError("Dimensions {} and {} of array inputs must match.".format(X.shape[1],box_size.shape[0]))

    # Check format of input
    Y = np.zeros(X.shape)
    for d in range(X.shape[1]):
        if np.any((X[:,d] > 3/2*box_size[d])|(X[:,d] <  - box_size[d]/2)):
            warnings.warn("Data positions must fit in box of dimension twice the one provided, otherwise the translation may lead to points out of the periodic box. One can wrap the data to fix this issue.")
        Y[:,d] = X[:,d] + box_size[d]/2
        Y[:,d][Y[:,d]>box_size[d]] -= box_size[d]
    return Y

def fuse_labels(l1,l2,debug=False):
    '''
    Find a union labelling systems between two independant labelling systems
     of the same size.

    TODO : Simplify/rewrite this function ... maybe using graph tools in python.
    '''
    # Check size of labels
    size = l1.shape[0]
    assert size == l2.shape[0], "Labels must be numpy arrays of the same size."
        
    # Check format of labels
    for labels in (l1,l2):
        assert np.min(labels) > -2 , "Incorrect format of labels, labels must be integer numbers, and outliers are set to -1."
        assert labels.dtype == int , "Incorrect format of labels, labels must be integer numbers, and outliers are set to -1."
    
    # Convert one label names into alphabet letters
    set1 = set(l1)
    set2 = set(l2)
    if set1.intersection(set2) is not set():
        ll2 = [chr(ord('`')+l+1) for l in l2]
    
    # Double labelling
    L = np.array([[label1,label2] for (label1,label2) in zip(l1,ll2)])
    
    # Initialize final set and labels
    S = {}
    all_labels = set()
    FL = np.zeros(size,dtype=int)
    
    # reduce size by excluding treatment of outliers
    L =   L[(l1!=-1) & (l2!=-1)]
    FL = FL[(l1!=-1) & (l2!=-1)]
    
    # Loop on each position
    for i,(label1,label2) in enumerate(L):
        # If the first label is already assigned in final set but not the other label
        if (label1 in all_labels) & (label2 not in all_labels):
            for k,labels_mix in S.items():
                # When the label is already assigned to one of the final label set
                if (label1 in labels_mix):
                    # add the correspondant label in the same cluster
                    labels_mix.add(label2)
                    all_labels.add(label2)
                    FL[i] = k
                    # No need to finish the loop on final set, since the label appears only once
                    break 
        # inverse case 
        elif (label2 in all_labels) & (label1 not in all_labels):
            for k,labels_mix in S.items():
                # When the label is already assigned to one of the final label set
                if (label2 in labels_mix):
                    # add the correspondant label in the same cluster
                    labels_mix.add(label1)
                    all_labels.add(label1)
                    FL[i] = k
                    # No need to finish the loop on final set, since the label appears only once
                    break
        # if both labels are already assigned in the final set
        elif (label2 in all_labels) & (label1 in all_labels):
            # look for the minimal index of final set to keep
            index = []
            for k,labels_mix in S.items():
                if (label1 in labels_mix)|(label2 in labels_mix):
                    index.append(k)
            # if two disctinct set are detected, copy all labels from the second set to the first one and delete the second one
            if len(index) != 1:
                S[index[0]] = S[index[0]].union(S[index[1]])
                S.pop(index[1],None)
                FL[i] = index[0]
                FL[FL == index[1]] = index[0]
            # if both labels are already assign in the same set, only update the label array
            else :
                FL[i] = index[0]
        
        # if none of the two labels are already assigned, create a new label in the final set
        else:
            k = 0 if len(S) == 0 else max(S.keys()) + 1
            FL[i] = k
            S[k] = {label1,label2}
            all_labels.add(label1)
            all_labels.add(label2)
            
    # Mapping of original labels to final ones
    map_start = {}
    for k,labels_mix in S.items():
        for el in labels_mix:
            map_start[str(el)] = k
    map_start['-1'] = -1
    map_start['`'] = -1

    # Correct labels of outliers
    FL_wo = np.ones(size,dtype=int)*-1
    FL_wo[(l1!=-1) & (l2!=-1)] = FL
    
    # outliers in the second labels
    for i,is_outlier_l2 in enumerate(l2 == -1):
        # If fisrt labels that are outliers in second labelling is not assigned, create a new final label
        if (is_outlier_l2 == False) | (str(l1[i]) in map_start.keys()): 
            continue
        k = 0 if len(S) == 0 else max(S.keys()) + 1
        FL_wo[i] = k
        S[k] = {str(l1[i])}
        all_labels.add(str(l1[i]))
        map_start[str(l1[i])] = k
        
    # outliers in the first labels
    for i,is_outlier_l1 in enumerate(l1 == -1):
        # If fisrt labels that are outliers in second labelling is not assigned, create a new final label
        if (is_outlier_l1 == False) | (str(ll2[i]) in map_start.keys()): 
            continue
        k = 0 if len(S) == 0 else max(S.keys()) + 1
        FL_wo[i] = k
        S[k] = {str(ll2[i])}
        all_labels.add(str(ll2[i]))
        map_start[str(ll2[i])] = k
    

    FL_wo[l2 == -1] = [map_start[str(el)] for el in  l1[l2 == -1]]
    FL_wo[l1 == -1] = [map_start[str(el)] for el in np.array(ll2)[l1 == -1]]

    # Check that all input labels in reduced ensemble (without outliers) have been assigned 
    set1.discard(-1)
    set2.discard(-1)
    assert np.array([len(S[j]) for j in S.keys()]).sum() == len(set1) + len(set2) , "All labels from inputs have not been assigned to an output label."
    
    # Map final labels to the smallest integer values
    map_min = {key:nl  for nl,key in enumerate(S.keys())}
    map_min[-1] = -1
    FL_wo  = np.array([ map_min[el] for el in FL_wo ])
    
    if debug:
        print("data first labelling :\n{}".format(l1),flush=True)
        print("data second labelling :\n{}".format(ll2),flush=True)
        print("data final labelling :\n{}".format(FL_wo),flush=True)
        print("Image Set of labels :\n{}".format(S),flush=True)
        print("First mapping all labels-> last:\n{}".format(map_start),flush=True)
        print("Second mapping in integer ensemble :\n{}".format(map_min),flush=True)
    return FL_wo
 
def analysis_cluster(ag,ag_split_by_atom,nmol,natoms,min_samples,eps,
                     with_distance_matrix=False,matrix_to_allocate=None,
                     with_double_labelling=False,debug=False,verbose=True):
    '''
    Get label clusters using different methods :
    
    1) with_distance_matrix = True and with_double_labelling = False
    Compute the distance matrix using periodic box condition usgin MDAnalysis.
    Pass a preallocated matrix array in matrix_to_allocate to reuse prexistent memory allocation.
    Slow method if the number of atoms overpass 5000. 
    
    2) with_distance_matrix = False and with_double_labelling = False
    Do not use periodic boundary conditions.
    Fast since using sklearn tools.
    
    3) with_distance_matrix = False and with_double_labelling = True (DEFAULT)
    First, the algorithm will first compute two labels using DBSCAN on sparse distance matrix 
    with radius neigbor pre-processing on both current positions and positions after half-diagonal box . 
    Then, a minimal labelling is set to each atom using these two labelling system.
    Only works with cubic box.
    Fast.
    
    TODO : implement other box shapes.
    
    '''
    centers = [mol.center_of_mass() for mol in ag_split_by_atom]
    center_array = np.array(centers)
    
    # Method 1
    if with_distance_matrix :
        if verbose: print("Method with total distance matrix",flush=True)
        if matrix_to_allocate is not None:
            assert matrix_to_allocate.shape[0] == int(nmol*natoms*(nmol*natoms - 1)/2) , 'Pre-allocated array must be an array with shape n*(n-1)/2'
        distance_matrix = self_distance_array(center_array, box=ag.dimensions,result=matrix_to_allocate)
        dbscan = DBSCAN(metric="precomputed",eps=eps,min_samples=min_samples)  
        all_cluster_labels = dbscan.fit_predict(squareform(distance_matrix))
        labels = all_cluster_labels.copy()
        labels_pbc = all_cluster_labels.copy()
        
    # Method 2
    elif not with_double_labelling:
        if verbose: print("No pbc considerations",flush=True)
        dbscan = DBSCAN(eps=eps,min_samples=min_samples)  
        all_cluster_labels = dbscan.fit(center_array).labels_
        labels = all_cluster_labels.copy()
        labels_pbc = all_cluster_labels.copy()
        
    # Method 3 
    else :
        if verbose: print("Method of double labelling",flush=True)
        # Compute sparse radius matrix distance with original positions 
        box = ag.universe.dimensions[:3]
        centers = [mol.center_of_mass() for mol in ag_split_by_atom]
        center_array = np.array(centers)
        matrix_distance = radius_neighbors_graph(center_array,radius=eps, mode='distance')

        # Compute sparse radius matrix distance on shifted positions 
        matrix_distance_pbc = radius_neighbors_graph(shift_pbc(center_array,box),radius=eps, mode='distance')
        
        # Clustering analysis
        db = DBSCAN(metric="precomputed",min_samples=min_samples,eps=eps).fit(matrix_distance)
        db_pbc = DBSCAN(metric="precomputed",min_samples=min_samples,eps=eps).fit(matrix_distance_pbc)

        # Find best labelling
        labels = db.labels_
        labels_pbc = db_pbc.labels_
        all_cluster_labels = fuse_labels(labels,labels_pbc,debug=debug)
        
        if verbose: print("# first labels :{},# second labels :{}".format(
            np.unique(labels).shape,np.unique(labels_pbc).shape),flush=True)

    (unique, counts) = np.unique(all_cluster_labels, return_counts=True)
    
    if verbose: print("# final labels{}".format(unique.shape),flush=True)
    
    # list of resids in the same order of count numbering
    ar_resids=[]
    for label in unique:
        ar_resids+=list(np.unique(ag.resids[all_cluster_labels==label]))
    
    if verbose:
        print("nb of resids labelled {} : nb of molecules : {}".format(len(ar_resids),nmol),flush=True)
        
    assert len(ar_resids)==nmol, 'Not every molecules are assigned to a unique cluster,\
    maybe it means that resids are not assigned properly.'

    # outliers are labeled -1
    if unique[0] == -1:
        # ignore outliers
        unique = unique[1:]
        counts = counts[1:]

    frequencies = np.asarray((unique, counts/natoms))
    frequencies_extend = np.zeros((3,nmol),dtype=np.int64)
    frequencies_extend[:frequencies.shape[0], :frequencies.shape[1]] = frequencies    
    frequencies_extend[2,:] = np.asarray(ar_resids)
    
    return frequencies_extend, labels , labels_pbc, all_cluster_labels

def get_infos(u,selection):
    '''
    Return number of molecules and number of atoms per molecule of solute.
    '''
    assert len(u.select_atoms(selection))>0 ,\
     "No atoms in the solute selection : {}".format(selection)

    nmol = np.unique(u.select_atoms(selection).resids).shape[0]
    natoms_l=[]
    for resid in u.select_atoms(selection).resids:
        natoms_l.append(u.select_atoms("{} and resid {}".format(selection,resid)).n_atoms)
    natoms = np.unique(natoms_l)
    assert len(natoms)==1 , "Molecules of the same type have not the same number of atoms ! Maybe redefine properly resid in the reference file. "
    natoms = natoms_l[0]
    return nmol,natoms

def lflat(lst):
    '''
    Flat a list of list into 1D list.
    '''
    lout=[]
    for el in lst :
        for eli in el:
            lout.append(eli)
    return lout

def get_resids_labels_from_cluster(sim,cluster_dict,frame):
    '''
    Return a numpy array of size (2,nmol) with labelled resids per aggregate cluster ID
    '''
    ind_resid=0
    ar_out=np.zeros((2,cluster_dict[sim].nmol),dtype=int)
    for labelID,count in cluster_dict[sim].counts[frame,0:2,:].T:
        #print(labelID,count,flush=True)
        start=ind_resid
        end=ind_resid+count
        ar_resids=cluster_dict[sim].counts[frame,2,start:end]
        ar_labelID=[labelID]*count
        ar_out[:,start:end]=np.array((ar_resids,ar_labelID))
        ind_resid+=count
    return ar_out

def save_input(args):
    logfile = os.path.basename(sys.argv[0]).split('.')[0]+'.log'
    with open('./{}/{}'.format(args.outdirname,logfile), 'w') as f:
        for key,value in vars(args).items():
            print(key,value,file=f)

def write_csv_properties(args,cluster_dict,sim_name,csvfile):
    '''
    Save number of cluster and maximal size of cluster aggregates in one CSV file.
    '''
    # Check consistency of data over several simulations
    sim0=sim_name[0]
    for sim in sim_name[1:]:
        assert cluster_dict[sim].times.shape==cluster_dict[sim0].times.shape,"All data must be recorded with the same number of frames."
    # Create the MultiIndex
    properties=['nclust','maxclust']
    midx = pd.MultiIndex.from_product([sim_name, properties])
    # Create sample data for each property, and add the MultiIndex.
    n_times=cluster_dict[sim0].times.shape[0]
    n_sim=len(sim_name)
    n_prop=len(properties)
    data_midx=np.array([[getattr(cluster_dict[sim],properties[i]) for i in range(n_prop)] for sim in cluster_dict]).reshape(2*n_sim,n_times)    
    df_prop = pd.DataFrame(data_midx, index = midx,columns=cluster_dict[sim0].times)
    df_prop.index.names = ['simulation','property']
    df_prop.columns.names=['time (ns)']    
    # Write to CSV
    df_prop.to_csv("{}/{}".format(args.outdirname,csvfile),float_format='%.3f')

def write_csv_histograms(args,cluster_dict,sim_name,csvfile):
    '''
    Save histograms of cluster size in one CSV file.
    '''
    # List all data from each simulation and concatenate into multiindex dataframe
    list_df=[]
    for sim in sim_name: 
        list_df.append(pd.DataFrame(getattr(cluster_dict[sim],'histo_clust').T,columns=("n","P(n)")))
    df_histo=pd.concat(list_df, keys=[sim for sim in sim_name],axis=1).T
    df_histo.index.set_names(["simulation","property"],inplace=True)
    # Write to CSV
    df_histo.to_csv("{}/{}".format(args.outdirname,csvfile),float_format='%.3f')

def write_csv_resids(args,cluster_dict,sim_name,csvfile):
    '''
    Save statistics and resids of all cluster aggregates at each frame in one CSV file.
    '''
    sim0=sim_name[0]
    df_resids=[]
    for frame in range(cluster_dict[sim0].n_frames):
        range_sim=[[sim for i in range(cluster_dict[sim].nmol)] for sim in sim_name]
        range_resids=[get_resids_labels_from_cluster(sim,cluster_dict,frame)[0].tolist() for sim in sim_name]
        arrays = [lflat(range_sim),lflat(range_resids)]
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=["simulation", "resid"])
        df_sim = pd.DataFrame(lflat([list(get_resids_labels_from_cluster(sim,cluster_dict,frame)[1]) for sim in sim_name]),
                      index=index,columns=["label"])
        df_resids.append(df_sim)
    df_resids=pd.concat(df_resids,axis=1)
    df_resids.columns=cluster_dict[sim0].times
    df_resids.columns.names=['time (ns)']
    # Write to CSV
    df_resids.to_csv("{}/{}".format(args.outdirname,csvfile),float_format='%.3f')

def __init__(self):
    return None

def main():
    # Read command-line arguments 
    args = parse_args(required=True)
    check_version()
    check_method(args)
    check_warnings(args)
    ldirs = check_dirs(args)
    create_dir(args)

    # Save the input arguments
    save_input(args) 
    
    # Loop on selected simulation directories
    cluster_dict = {}
    for sim in ldirs:
        u = load_traj(sim,args)
        if args.selection:
            selection = args.selection
            nmol,natoms = get_infos(u,selection)
        else:
            selection = "resname {}".format(args.resname)
            nmol,natoms = get_infos(u,selection)
        ag = u.select_atoms(selection)
        cluster_dict[sim]=Clustsize(ag, verbose=False,nmol=nmol,natoms=natoms,
                                min_samples=args.mins,eps=args.eps,
                                with_distance_matrix=args.with_distance_matrix,
                                with_double_labelling=args.with_double_labelling
                                ).run(start=args.begin,stop=args.end,verbose=(args.verbose>1))

    # Write results
    write_csv_properties(args,cluster_dict,ldirs,csvfile='cluster_properties.csv')
    write_csv_histograms(args,cluster_dict,ldirs,csvfile='cluster_histo.csv')
    write_csv_resids(args,cluster_dict,ldirs,    csvfile='cluster_resids.csv')

if __name__ == "__main__":
    main()