# multiClustsize

Find all aggregate clusters for a set of MD trajectories.

Quick start
===========

Parameters
-------
``` bash
usage: run_Clustsize.py [-h] -i INPUTDIR -m RESNAME [-o OUTDIRNAME]
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

Outputs
-------
The **_run_Clustsize.py_** script generates three CSV files located in `OUTDIRNAME` directory (by default : `output`) and these data are easily editables using, for instance, the DataFrame format from the **pandas** library.
- **cluster_histo.csv** : size indices in columns and $simulation$ key, aggregation number $n$ and size distribution $P(n)$ in rows.
- **cluster_properties** : time (in ns) in columns and $simulation$ key, number of clusters $nclust$ and number of molecules in the biggest cluster $maxclust$ in rows.
- **cluster_resids.csv** : time (in ns) in columns and $simulation$ key and $resid$ indices in rows. For a given simulation, time and resid index one have stored the label (positive or ) of the cluster where is found the molecule. 

More insights of the output data is accesible in the test directory `./tests/test_analysis/`.


To do
-------
How to read easily these csv files using **pandas**.
Documentation for the use of the plotting script **_read_csv_and_plot.py_**.
