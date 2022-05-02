# Phenix with SPaDES demo
Demo to replicate the method from the submitted article on combining Phenix with SPaDES for improved structure refinement.

###

## REQUIREMENTS

In brief, Phenix and Rosetta with SPaDES must be installed before running the combined protocol. The approach described in our manuscript was tested with Rosetta 3.7 and Phenix version 1.19.2-4158

Phenix with SPaDES requires Python 2.7 to work. This can be done using conda environments, e.g.:

```
conda create --name py27 python=2.7
```

Creates a python 2.7 environment within GNU/Linux and other unix-based operating systems. This can then be activated:

```
conda activate py27
```

To run the demo. You can deactivate the environment with:

```
conda deactivate
```

And subsequently delete the environment if you desire:

```
conda remove --name py27 --all
```

Additional requirements (in general for combining phenix with Rosetta), include the following:

* Python2 dev tools


These can be installed with:
```
sudo apt-get install python-dev
```
(Or equivalent on other distributions)

## INSTALLATION

The test demo includes some small hard coded statements in Phenix to facilitate the inclusion of waters in refinement. In Rosetta, some hydration related movers have been added to ensure compatibility. These will both be merged into the main branches of Phenix and Rosetta soon, but for ease of testing we include our modified versions of Phenix and Rosetta here.

First you will need to get the tarballs from Zenodo, and move them into the phenix_with_spades_demo folder: https://zenodo.org/record/5784150#.YboWkCYo_Jk

Then, move them into this directory (phenix_with_spades_demo) and decompress them:

```
tar -xvf rosetta.tar.gz
tar -xvf phenix.tar.gz
```
Next, install Phenix:

```
cd phenix/
./install --prefix /path/to/phenix/phenix/install_dir
```

Where `/path/to/phenix/` is the current path to Phenix (e.g. `/home/lucas/phenix_with_spades_demo/phenix`)
After Phenix has installed, you will need to prepare your paths to compile the interfaced version of Rosetta with Phenix. Make sure you put the correct path locations in:

```
source /path/to/phenix/phenix/install_dir/phenix-1.19.2-4158/phenix_env.sh
ln -s /path/to/phenix/phenix/install_dir/phenix-1.19.2-4158/conda_base/lib/libpython2.7.so /path/to/rosetta/rosetta/main/source/external/lib/
LD_LIBRARY_PATH=/path/to/rosetta/rosetta/main/source/build/src/release/linux/5.13/64/x86/gcc/9/python
```

Note that you won't yet have actually built the rosetta/main/source/build/ folder, but without Rosetta knowing which folder to search for specific python libraries, it will fail on install. The `9` in the LD_LIBRARY_PATH location is depending on your version of gcc. One solution, if you're unsure of what to put, is to run the Rosetta install (as follows) without setting LD_LIBRARY_PATH. It will construct the folder, but fail. Then you can check the exact path needed, and set the LD_LIBRARY_PATH as appropiate. You can install Rosetta with:

```
cd rosetta/main/source
./scons.py bin mode=release -j 10 extras=python
```

Where `j` is the number of cores you want to install across at the same time. Note that this version of Rosetta will not search for GitHub version files by default (which is opposite to the nominal behaviour of Rosetta).

After Rosetta has installed, you will need to ensure you have the correct paths present to run the interfaced Rosetta:

```
export PHENIX_ROSETTA_PATH=/path/to/rosetta/rosetta/
```

Now you can finish interfacing Phenix with Rosetta:

```
rosetta.build_phenix_interface nproc=10
```

## RUNNING THE DEMO

Two folders are included that are critical to running the demo: `input_files` and `run_files`. You shouldn't need to edit any files in any of these locations. All the demo requires is listed below.

Update your environment variables further for the demo (if not following on from install):

```
source /path/to/phenix/phenix/install_dir/phenix-1.19.2-4158/phenix_env.sh
export PHENIX_ROSETTA_PATH=/path/to/rosetta/rosetta/
LD_LIBRARY_PATH=/path/to/rosetta/rosetta18/main/source/build/src/release/linux/5.8/64/x86/gcc/9/python
```

First, change into the run_files folder (this is where we will run the demo)

```
cd run_files
```

Next run prepare_input.sh in run_files (helps to prepare input files)

```
./prepare_input.sh
```

Note that if you move the Rosetta and Phenix folders around from their current relative locations, that can ruin the demo setup paths.

Now run the Phenix refine with SPaDES simulation with:

```
./refine.sh
```

Note this is setup to run on one processor. You can increase this with `nproc=2` (for 2 processors), but this is a rather memory heavy calculation, so be careful how many processors you run on, or run on a cluster. This calculation (on one processor) can take up to 6 hours.

At the end of the calculation, you should have a new folder called rosetta_1 (or 2, 3 etc. if you have run multiple times). Inside this folder, you will find the output of phenix.refine, and a folder called hydrate_structure. Inside this folder, you will find lilly_refined_final.pdb, which is the final structure of interest post refinement.
