# multiRDF

Run sequentially Radial Distribution Function analysis from MD trajectories.

Quick start
===========


Installation
------------
Here is the steps to install multiRDF in `~/programs`, but you can of course change the installation directory.

1. Clone this repository on your computer.
   ```bash
   mkdir -p ~/programs
   cd ~/programs
   git clone https://gitlab.galaxy.ibpc.fr/hardiagon/multiRDF.git
   ```

2. Add main script `run_RDF.py` to your python modules. 
For example, if you have a directory with some executable python script like /path/to/python/modules/, you need to link the main script on this program into this directory:
    ```bash
   ln -s ~/programs/multiRDF/run_RDF.py /path/to/python/modules/ 
   ```
3. Install the dependendies for Python.
For example, with conda, tou can create a new environment with :
    ```bash
    conda env create --file "~/programs/multiRDF/multirdf.yml" -n multirdf
    conda activate multirdf
   ```

Tutorial
--------

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
Then you need to specify the parameters you want to change for each analysis :

    TASK  SIM     OUTDIRNAME  FREQ   SEL1                       EXCLUSION_BLOCK   
    taks1 U01_3ch run0        1000  "resname U01 and name N12"  (1,1)
    task2 U02_3ch run0        1000  "resname U02 and name N4"   (1,1)

Then you just have to run the main script in a directory where you want to store results:
```
run_RDF.py
```

Parameters
-------
TODO : explanation on all parameters 

Outputs
-------
TODO : describe architecture of the output directories and output file formats