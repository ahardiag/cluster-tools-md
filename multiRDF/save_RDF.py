#!/usr/bin/env python 

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF,InterRDF_s
# https://docs.mdanalysis.org/1.0.0//documentation_pages/analysis/rdf.html
import numpy as np
import os, sys, time, glob, re ,math
import argparse
import pandas as pd
from argparse import RawTextHelpFormatter
import matplotlib.pyplot as plt

def parse_args(required=False):
    parser = argparse.ArgumentParser(description="Compute Radial Distribution Function using MDAnalysis module\n"
                                            "For selection commands, see MDAnalysis doc :\n"
                                            "https://docs.mdanalysis.org/stable/documentation_pages/selections.html",
                                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--traj",
                        help="Trajectory file (.xtc,.trr).")
    parser.add_argument("-r", "--ref",
                        required=required,
                        help="Reference or first frame file (.pdb, .gro).")
    parser.add_argument("--sel1",
                        required=required,
                        help="AtomGroup selection text.")
    parser.add_argument("--sel2",
                        help="AtomGroup selection text. Default: sel1 value")
    parser.add_argument("-o", "--outdirname",
                        default=".",
                        help="Output Sub-directory name.")
    parser.add_argument("-x", "--suffixname",
                        help="Output suffix name.")
    parser.add_argument("-b", "--begin",
                        default=0,
                        type=int,
                        help="Start frame of analysis.")
    parser.add_argument("-e", "--end",
                        type=int,
                        help="Stop frame of analysis.")
    parser.add_argument("--sub",
                        default=1,
                        type=int,
                        help="Subdivision of the trajectory, number of frames per window.\n"
                        "This will generated multiple traces on the same graph.")
    parser.add_argument("-m", "--max",
                        type=int,
                        default=15,
                        help="Maximum radial distance in Angstrom.")
    parser.add_argument("--ymax",
                        type=int,
                        default=20,
                        help="Maximum Y value.")
    parser.add_argument("--nbins",
                        type=int,
                        default=75,
                        help="Number of points of the output data.")
    parser.add_argument("--exclusion_block",
                        default=None,
                        nargs="+",
                        type=int,
                        help="Tuple with the number of atoms per molecule.\n"
                        "It avoids to count atoms in the same molecule.")
    parser.add_argument("-c","--cumulative",
                        action='store_true',
                        help="Returns cumulative number RDF, i.e. the number\n"
                        " i.e. the average number of particles within a distance r.\n"
                        "Only has a meaning if exclusion_block=(1,1)")
    parser.add_argument("-v","--verbose",
                        default=True,
                        help="Turn on progress bar for single RDF calculation.")
    parser.add_argument("--fontsize",
                        default=10,
                        help="Fontsize for legend, ticks and axis labels.")

    #args = parser.parse_args()
    return parser.parse_args()

def run(args):
    if args.exclusion_block is not None:
        ltuple=len(args.exclusion_block)
        assert ltuple < 3 , 'Error in input --exclusion_block, must be one or two integers !'
        if ltuple==1:
            args.exclusion_block=(args.exclusion_block,args.exclusion_block)
        elif ltuple==2:
            args.exclusion_block=tuple(args.exclusion_block)

    if args.sel2 is None:
        args.sel2 = args.sel1

    # Create output directory if not exist
    if args.verbose:
        try:
            os.makedirs("Figures")
            print("Directory ./Figures Created ")
        except FileExistsError:
            pass
        #    print("Directory ./Figures already exists")  
        try:
            os.makedirs("Output")    
            print("Directory ./Output Created ")
        except FileExistsError:
            pass
        #   print("Directory ./Output already exists") 

    # Check existence of trajectory files and load them
    assert os.path.exists(args.ref), 'File {} do not exist !'.format(args.ref)
    
    if args.traj is None:
        u=mda.Universe(args.ref)
    else:
        assert os.path.exists(args.traj), 'File {} do not exist !'.format(args.traj)
        # Load trajectory
        u=mda.Universe(args.ref,args.traj)

    # Check end time for analysis is smaller than simulation length time

    if args.end is None :
        args.end = (u.trajectory.n_frames*u.trajectory.dt/1000)
    elif (args.end*1000/u.trajectory.dt > u.trajectory.n_frames):
        args.end = (u.trajectory.n_frames*u.trajectory.dt/1000)

    # Output files
    if args.suffixname is None:
        args.suffixname="%s_%s_%d_%d"%(args.sel1.replace(" ","_"),args.sel2.replace(" ","_"),args.begin,args.end)

    try:
        os.makedirs("./Output/{}".format(args.outdirname))
        os.makedirs("./Figures/{}".format(args.outdirname))
    except:
        pass

    if args.traj is None:    
        filename=os.path.basename(args.ref).split('.')[0]
    else:
        filename=os.path.basename(args.traj).split('.')[0]
    
    outcsv = './Output/{}/rdf_{}.csv'.format(args.outdirname,args.suffixname)
    outpng ='./Figures/{}/rdf_{}.png'.format(args.outdirname,args.suffixname)

    # AtomGroup selections
    ag1=u.select_atoms(args.sel1)
    ag2=u.select_atoms(args.sel2)
    assert (ag1.n_atoms)*(ag2.n_atoms) > 0, 'An atom selection is empty, correct input --sel1 or/and --sel2'

    # Compute RDF
    ltime=[int(args.begin+k*(args.end-args.begin)/args.sub) for k in range(args.sub+1)]
    #print(ltime) #DEBUG    
    fig,ax=plt.subplots()
    df=pd.DataFrame(np.linspace(0,args.max,args.nbins),columns=['r_A'])

    for (itime,etime) in zip(ltime[:-1],ltime[1:]):

        # convert time if necessary
        iframe = int(itime*1000/u.trajectory.dt)
        if etime==0:
            # one frame case
            eframe=1
        else:
            eframe = int(etime*1000/u.trajectory.dt)
        
        # print(iframe,eframe) #DEBUG
        rdf = InterRDF(ag1,ag2,nbins=args.nbins,range=(0,args.max),
                       exclusion_block=args.exclusion_block,
                       verbose=args.verbose)
        rdf.run(start=iframe,stop=eframe,verbose=True)

        # Save data in CSV format (Output directory)
        # and returns cumulative number RDF if option `-c` requested
        cname="RDF_{}_{}".format(itime,etime)
        if args.cumulative:        
            cdf = np.cumsum(rdf.results.count)/(eframe-iframe)/ag1.n_atoms
            df[cname] = cdf    
        else:
            df[cname]=rdf.results.rdf
        # Plot graphs
        ax.plot(rdf.results.bins[1:], df[cname][1:],
                '-',label='[%d,%d] ns'%(itime,etime))
        #ax.plot(rdf.results.bins[1:], rdf.results.rdf[1:],
        #        '-',label='[%d,%d] ns'%(itime,etime))

    df.to_csv(outcsv,sep=' ',index=False)
    plt.xlim(0,args.max)
    plt.ylim(0,args.ymax)
    plt.xlabel(r'r($\AA$)',fontsize=args.fontsize)
    plt.ylabel('g(r)',fontsize=args.fontsize)
    plt.legend(ncol=2,fontsize=args.fontsize)
    plt.xticks(fontsize=args.fontsize)
    plt.yticks(fontsize=args.fontsize)
    fig.tight_layout()
    #fig.show() #DEBUG
    plt.savefig(outpng,dpi=100)

def __init__(self):
    return None

def main():
    args=parse_args(required=True)
    run(args)

if __name__ == "__main__":
    main()

