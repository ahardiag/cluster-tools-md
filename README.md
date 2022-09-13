# cluster-tools-md

Command-line analysis routines of Molecular Dynamics data with aggregating molecules.

## Installation

Here are the steps to install **cluster_tools** in `~/programs`, but you can of course change the installation directory.

1. Clone this repository on your computer.
   ```bash
   mkdir -p ~/programs
   cd ~/programs
   git clone https://github.com/ahardiag/cluster_tools_md
   ```

2. Use a python environment to install the package on it, e.g. using `conda` :
   ```bash
    conda create -n cluster_tools_md python=3.10
   ```

3. Build the package using `pip` :
   ```bash
   cd ~/programs
   pip install ./cluster_tools_md
   ```

### Usage
```bash
clustsize -h
eccentricity -h
multirdf -h
radialdens -h
```
### Development

To retrieve the same conda environment used during the development, use :
   ```bash
   conda env create --file "~/programs/cluster_tools/cluster_tools_md.yml" -n cluster_tools_md
   conda activate cluster_tools_md
   ```


-------

## clustsize

Computes the **size distributions** of atomic clusters in gas phase for a set of MD trajectories.
This is equivalent to `gmx clustsize` in GROMACS, but provides more options and editable outputs.

### Parameters

``` bash
usage: clustsize [-h]    -i INPUTDIR -m RESNAME [-o OUTDIRNAME]
                        [-r REFFILE] [-f TRJFILE] [-b BEGIN] -e END [-v]
                        --mins MINS --eps EPS [--with-distance-matrix]
                        [--with-double-labelling] [--ignore-warnings]
                        [--pattern-dirs PATTERN_DIRS]

Find all aggregate clusters for a set of MD trajectories.
For selection commands, see MDAnalysis doc :
https://docs.mdanalysis.org/stable/documentation_pagses/selections.html

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTDIR, --inputdir INPUTDIR
                        directory with all trajectories to analyse
  -m RESNAME, --resname RESNAME
                        Resname or molecule name. e.g: HC8
  -o OUTDIRNAME, --outdirname OUTDIRNAME
                        Output Sub-directory name.
  -r REFFILE, --reffile REFFILE
                        Commmon name of reference file in all simulation directories.
  -f TRJFILE, --trjfile TRJFILE
                        Commmon name of trajectory file in all simulation directories.
  -b BEGIN, --begin BEGIN
                        Start frame of analysis.
  -e END, --end END     Stop frame of analysis.
  -v, --verbose         Turn on for more details during the calculation.
  --mins MINS           Minimum number of points (atoms) to define a cluster.
  --eps EPS             Maximum distance between two points (atoms) in the same cluster.
  --with-distance-matrix
                        Use total matrix distance to treat periodic boundary conditions. Slow for large systems!
  --with-double-labelling
                        Use double labelling method to perform clustering analysis on two shifted set of positions. Allow to search for clusters close to the box boundaries.
  --ignore-warnings     Avoid stdout all warnings from using clustering on sparse matrix.
  --pattern-dirs PATTERN_DIRS
                        String pattern to select a serie of data directories with similar pattern in the input directory.
```

### Outputs

The script generates three CSV files located in `OUTDIRNAME` directory (by default : `output`) and these data are easily editables using, for instance, the DataFrame format from the **pandas** library.
- **cluster_histo.csv** : size indices in columns and $simulation$ key, aggregation number $n$ and size distribution $P(n)$ in rows.
- **cluster_properties** : time (in ns) in columns and $simulation$ key, number of clusters $nclust$ and number of molecules in the biggest cluster $maxclust$ in rows.
- **cluster_resids.csv** : time (in ns) in columns and $simulation$ key and $resid$ indices in rows. For a given simulation, time and resid index one have stored the label (positive or ) of the cluster where is found the molecule. 

More insights of the output data is accesible in `tests` directory.

-------

## multirdf

Run sequentially Radial Distribution Function analysis for a set of MD trajectories.

### Tutorial

First you need to create a config file with some default parameters, in a file on the root location : `~/.runRDFrc`:
```ini
  [DefaultParameters]
    BEGIN=              100             # Start time in ns for the analysis
    END=                800
    MAX=                20              # Maximum radial distance (in
                                        # Angstrom)
    SUB=                2               # Number of time intervals to
                                        # divide the
                                        # Trajectory total time 
    NBINS=              200             # Number of points for one RDF 
                                        # trace        
    PATH=               ../Data_traj/   # Relative path to data from the
                                        # current directory 
    EXCLUSION_BLOCK=    (1,1)           
```
Then you need to specify the parameters you want to change for each analysis in a input file, e.g. `parameters.in`:

    TASK  SIM     OUTDIRNAME  FREQ   SEL1                       EXCLUSION_BLOCK   
    taks1 U01_3ch run0        1000  "resname U01 and name N12"  (1,1)
    task2 U02_3ch run0        1000  "resname U02 and name N4"   (1,1)

Then you just have to run the main script in a directory where you want to store results:
```bash
multirdf
```

--------

## eccentricity

Compute the eccentricity of each cluster found after using `clustsize` command.

```bash
usage: eccentricity [-h] -s SIM -i INPUTDIR -c CSVDIR [--nbins NBINS] [-m RESNAME] [-o OUTDIRNAME]
                    [--outfigdirname OUTFIGDIRNAME] [-r FILEREF] [-f FILETRJ]
                    [--method {double_loop,best_centering}] [--compound {residues,atoms}]
                    [--rename RENAME] [--in_memory] [-z SIZE [SIZE ...]] [--range] [-v]

Compute Spherical Radial Density around cluster COM.

options:
  -h, --help            show this help message and exit
  -s SIM, --sim SIM     Name of the Data trajectory.
  -i INPUTDIR, --inputdir INPUTDIR
                        directory with all trajectories to analyse.
  -c CSVDIR, --csvdir CSVDIR
                        directory with the dataframe `cluster_resids.csv`.
  --nbins NBINS         Number of bins of the density distribution.
  -m RESNAME, --resname RESNAME
                        Resname or molecule name. e.g: HC8.
  -o OUTDIRNAME, --outdirname OUTDIRNAME
                        Output Sub-directory name.
  --outfigdirname OUTFIGDIRNAME
                        Output Sub-directory name for figures.
  -r FILEREF, --fileref FILEREF
                        Commmon name of reference file in all 
                        simulation directories.
  -f FILETRJ, --filetrj FILETRJ
                        Commmon name of trajectory file in all 
                        simulation directories.
  --method {double_loop,best_centering}
                        Choose the method for computing the inertia tensor.
  --compound {residues,atoms}
                        Choose the level at which is computed the inertia tensor.
                        `residues` is faster as it only consider the center of mass of the molecules.
  --rename RENAME       Set to True to consider sim name with no suffix `Data_`.
  --in_memory           Charge the whole trajectory in RAM memory.
                        Faster than reading from file.
  -z SIZE [SIZE ...], --size SIZE [SIZE ...]
                        Choose specific size(s) of clusters.
  --range               Range of values provided by three values `start end step`by -z/--size option.
  -v, --verbose         Turn on for more details during the 
                        calculation.
```

--------

## radialdens

Compute the radial density of an atom group around each cluster, depending on its size.

```bash
usage: radialdens [-h] -s SIM -z SIZE [SIZE ...] [--range] -i INPUTDIR -c CSVDIR [--nbins NBINS]
                  [-m RESNAME] [-o OUTDIRNAME] [--outfigdirname OUTFIGDIRNAME] [-r FILEREF]
                  [-f FILETRJ] [--rename RENAME] [--in_memory] [-v]

Compute Spherical Radial Density around cluster COM.

options:
  -h, --help            show this help message and exit
  -s SIM, --sim SIM     Name of the Data trajectory.
  -z SIZE [SIZE ...], --size SIZE [SIZE ...]
                        Size of the cluster.
  --range               Range of values provided by three values `start end step`by -z/--size option.
  -i INPUTDIR, --inputdir INPUTDIR
                        directory with all trajectories to analyse.
  -c CSVDIR, --csvdir CSVDIR
                        directory with the dataframe `cluster_resids.csv`.
  --nbins NBINS         Number of bins of the density distribution.
  -m RESNAME, --resname RESNAME
                        Resname or molecule name. e.g: HC8.
  -o OUTDIRNAME, --outdirname OUTDIRNAME
                        Output Sub-directory name.
  --outfigdirname OUTFIGDIRNAME
                        Output Sub-directory name for figures.
  -r FILEREF, --fileref FILEREF
                        Commmon name of reference file in all 
                        simulation directories.
  -f FILETRJ, --filetrj FILETRJ
                        Commmon name of trajectory file in all 
                        simulation directories.
  --rename RENAME       Set to True to consider sim name with no suffix `Data_`.
  --in_memory           Charge the whole trajectory in RAM memory.
                        Faster than reading from file.
  -v, --verbose         Turn on for more details during the 
                        calculation.
```

--------

## Testing

A set of scripts (`test1.sh`,`test2.sh`,...) are given as examples and control tests, it processes the data from the article cited in section [Citation](#Citation).

### Usage
Run the two first tests
```bash
cd ~/programs/cluster_tools_md/tests/
./run_tests.sh test1 test2
```

--------

## Citing

If you use **cluster_tools** in your research, please cite the following article:

[coming soon]

