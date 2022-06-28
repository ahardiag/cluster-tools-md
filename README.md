# cluster_tools

Python scripts to analyze Molecular Dynamics simulations with a phase that forms clusters.

Quick start
===========


Installation
------------
Here are the steps to install **cluster_tools** in `~/programs`, but you can of course change the installation directory.

1. Clone this repository on your computer.
   ```bash
   mkdir -p ~/programs
   cd ~/programs
   git clone https://gitlab.galaxy.ibpc.fr/hardiagon/cluster_tools.git
   ```

2. Add all python scripts to your python modules.
For example, if you have a directory with some executable python script like /path/to/python/modules/, you need to link the main script on this program into this directory:
    ```bash
   ln -s ~/programs/cluster_tools/multiclustsize/multiclustsize.py /path/to/python/modules/ 
   ```
3. Install the dependencies :
These programs use the python packages `mdanalysis`, `pandas`, `numpy`,`numba`, or `scikit-learn` :

For example, with conda, you can create a new environment with all the python packages, using the one provided in this directory :
   ```bash
   conda env create --file "~/programs/cluster_tools/cluster_tools.yml" -n cluster_tools
   conda activate cluster_tools
   ```


What for ?
--------

- `multiclustersize/run_Clustsize.py` computes the size distributions of atomic clusters in gas phase.
This is equivalent to `gmx clustsize` in GROMACS, but provides more options and editable outputs.
- `multiRDF` is a wrapper of MDAnalysis tools for computing Radial Distribution Functions from MD trajectories.
In particular, it allows to apply sequentially the same analysis on several trajectories.

See the README file in all subdirectories.


Citation
-------
If you use **cluster_tools** in your research, please cite the following article:


