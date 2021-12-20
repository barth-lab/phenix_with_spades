# Phenix with spades demo
Demo to replicate results from submitted article on combining phenix with spades for improved structure refinement

###

## REQUIREMENTS

In brief, Phenix and Rosetta with SPaDES must be installed before running the combined protocol. The approach described in our manuscript was tested with Rosetta 3.7 and Phenix version 1.19.2-4158

Phenix with spades requires Python2.7 to work. This can be done using conda environments, e.g.:

```conda create --name py27 python=2.7```

Creates a python 2.7 environment within linux. This can then be activated:

```conda activate py27```

To run the demo. You can deactivate the environment with:

```conda deactivate```

And subsequently delete the environment if you desire:

```conda remove --name py27 --all```

## INSTALLATION

The test demo includes some small hard coded statements in phenix to facilitate the inclusion of waters in refinement. In Rosetta, some hydration related movers have been added to ensure compatibility. These will both be merged into the main branches of phenix and rosetta soon, but for ease of testing we include our modified versions of phenix and Rosetta here.

First you will need to get the tarballs from zenodo, and move them into the phenix_with_spades_demo folder: https://zenodo.org/record/5784150#.YboWkCYo_Jk

Then, move them into this directory (phenix_with_spades_demo) and decompress them:

```tar -xvf rosetta.tar.gz
tar -xvf phenix.tar.gz```

Next, install rosetta:

```cd rosetta/main/source
./scons.py bin mode=release -j 10 extras=python```

Where `j` is the number of cores you want to install across at the same time.

Then, install phenix with:

```cd phenix/
./install --prefix /path/to/phenix/phenix/install_dir```

Where `/path/to/phenix/` is the current path to phenix (e.g. `/home/lucas/phenix_with_spades_demo/phenix`)

After phenix has installed, you need to compile the interfaced version of rosetta with phenix:

```source /path/to/phenix/phenix/install_dir/phenix-1.19.2-4158/phenix_env.sh
export PHENIX_ROSETTA_PATH=/path/to/rosetta/rosetta/
LD_LIBRARY_PATH=/path/to/rosetta/rosetta18_UCS4/main/source/build/src/release/linux/5.8/64/x86/gcc/9/python```

This last line is to ensure there is no confusion in case you have multiple versions of rosetta installed. Note some of the numbers in this line (e.g. 9 towards the end), refers to your gcc compiler, so may be slightly different.

Now you can install phenix with rosetta:

```rosetta.build_phenix_interface nproc=10```

Note, you may get an error immediately stating that:

> phenix.python2.7 couldn't be found

To fix this, you must change:

```phenix.python2.7 options.py```

to 

```phenix.python options.py``` 

in update_options.sh, and:

```phenix.python2.7 update_ResidueType_enum_files.py```

to

```phenix.python update_ResidueType_enum_files.py```

in update_ResidueType_enum_files.sh

## RUNNING THE DEMO

Two folders are included that are critical to running the demo: `input_files` and `run_files`. You shouldn't need to edit any files in any of these locations. All the demo requires is listed below.

Update your environment variables further for the demo (if not following on from install):

```source /path/to/phenix/phenix/install_dir/phenix-1.19.2-4158/phenix_env.sh
export PHENIX_ROSETTA_PATH=/path/to/rosetta/rosetta/
LD_LIBRARY_PATH=/path/to/rosetta/rosetta18/main/source/build/src/release/linux/5.8/64/x86/gcc/9/python```

First, change into the run_files folder (this is where we will run the demo)

```cd run_files```

Next run prepare_input.sh in run_files (helps to prepare input files)

```./prepare_input.sh```

Note that if you move the rosetta and phenix folders around from their current relative locations, that can ruin the demo setup paths.

Now run the phenix refine with spades simulation with:

```./refine.sh```

Note this is setup to run on one processor. You can increase this with `nproc=2` (for 2 processors), but this is a rather memory heavy calculation, so be careful how many processors you run on, or run on a cluster. This calculation (on one processor) can take up to 6 hours.

At the end of the calculation, you should have a new folder called rosetta_1 (or 2, 3 etc. if you have run multiple times). Inside this folder, you will find the output of phenix.refine, and a folder called hydrate_structure. Inside this folder, you will find lilly_refined_final.pdb, which is the final structure of interest post refinement.
