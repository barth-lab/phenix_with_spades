# Phenix with SPaDES

This repo contains information pertaining to interfacing and running Phenix with Rosetta SPaDES for improved structure refinement based on our recent paper TBC.

We include specific installation instructions as well as detailed demos to demonstrate refinement.

In running the ion sampling and SPaDES designs, one only needs the most recent build of Rosetta (https://github.com/RosettaCommons/main) and the new hydrate functions provided in our corresponding repo (https://github.com/barth-lab/SPaDES_DESIGNS). 

We provide here **ALL** functions and movers needed to run the new SPaDES, including those specific functions needed for design/ion sampling given in the other repo, and the movers/protocols needed for phenix. These new movers can be used in any Rosetta Scripts XML context. These will soon be moved into the main branch of Rosetta.

## Requirements

In brief, Phenix and Rosetta with SPaDES must be installed before running the combined protocol. The approach described in our manuscript was initially tested with Rosetta 3.7 (2018) and Phenix version 1.19.2-4158, however we have since benchmarked the installation/refinement with the most up-to-date versions of Rosetta (weekly release rosetta.351, 17/07/2023) and Phenix (Phenix 1.20.1-4487) and obtained the same outcome.

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

Please note that, while it is not yet released, it would appear the Phenix developers are planning a Python3-based build since Python2 had its last major update in 2010. If you install this build there is no guarantee it will work with Rosetta SPaDES, though we will try to update when the time comes.

## Installation

We have added some additional flags passable to Phenix to facilitate the inclusion of waters in refinement as well as various new/updated hydrate Rosetta functions and movers that can be installed with the latest version of Rosetta. These will all be soon merged into the latest branches of both software. 

First, you will need to download the latest Phenix build (https://phenix-online.org/download) and Rosetta (https://www.rosettacommons.org/software/license-and-download).

Then, move them into a directory as you prefer and decompress them:
```
tar -xvf rosetta.tar.gz
tar -xvf phenix.tar.gz
```
Before installing phenix, you need to copy over the changes to the rosetta refinement files such that the spades arguments can be read in
```
cd /path/to/this/folder/phenix_changes/
cp refine.py /path/to/phenix/modules/phenix/phenix/rosetta/
cp rosetta_refine.py /path/to/phenix/modules/phenix/phenix/command_line/
```
Where `/path/to/phenix/` is the current path to Phenix (e.g. `/home/lucas/phenix`). When running ```phenix.rosetta_refine``` later on, one can access these extra flags by parsing the ```spades_args=/path/to/flags.file``` argument

Next we can install phenix
```
cd phenix/
./install --prefix /path/to/phenix/phenix/install_dir
```
After Phenix has installed, you will need to prepare your paths to compile the interfaced version of Rosetta with Phenix. Make sure you put the correct path locations in:
```
source /path/to/phenix/phenix/install_dir/phenix-1.20.1-4487/phenix_env.sh
export PHENIX_ROSETTA_PATH=/path/to/rosetta/
```
Prior to installing Rosetta, you will need to copy over the necessary movers/options/protocol files given here.
```
cd /path/to/this/folder/hydrate_movers
cp Hydrate* hydrate* Remove* /path/to/rosetta/main/source/src/protocols/hydrate/
cd OPTIONS/
cp options_rosetta.py /path/to/rosetta/main/source/src/basic/options/
cp site.settings /path/to/rosetta/main/source/tools/build/
cd ../PhenixInterface 
cp PhenixInterface.cc /path/to/rosetta/main/source/src/core/scoring/cryst/
cd ../PROTOCOL_FILES
cp protocols*.settings /path/to/rosetta/main/source/src/
cp init.Mover*ihh /path/to/rosetta/main/source/src/protocols/init/
```
While the hydrate, options, protocol and Mover.ihh files are introduced by us for SPaDES, the site.settings is a necessary insert from phenix (https://phenix-online.org/documentation/extras/site.settings) while the PhenixInterface correction was proposed by a user owing to missing imports (https://rosettacommons.org/node/11491)

Now you've copied the necessary files, we can install Rosetta:
```
cd rosetta/main/source
rosetta.build_phenix_interface nproc=10
```
Where nproc represents the number of cores you want to install across at the same time. The above command replaces the more standard
```
./scons.py bin mode=release -j 10
```

## Running the demo

Two folders are included that are critical to running the demo: `input_files` and `run_files`. You will need to make two small changes in `refine.sh` to indicate the path to Rosetta and its database (see `PATH/TO/ROSETTA` in file). After that, you shouldn't need to edit any files in any of these locations

Update your environment variables further for the demo (if not following on from install):
```
source /path/to/phenix/phenix/install_dir/phenix-1.19.2-4158/phenix_env.sh
export PHENIX_ROSETTA_PATH=/path/to/rosetta/rosetta/
```
Change into the run_files folder (this is where we will run the demo)
```
cd run_files
```
Next run prepare_input.sh in run_files (helps to prepare input files)
```
./prepare_input.sh
```
Now run the Phenix refine with SPaDES simulation with:
```
./refine.sh
```

Note this is setup to run on one processor. You can increase this with `nproc=2` (for 2 processors), but this is a rather memory heavy calculation, so be careful how many processors you run on, or run on a cluster. This calculation (on one processor) can take up to 6 hours.

At the end of the calculation, you should have a new folder called rosetta_1 (or 2, 3 etc. if you have run multiple times). Inside this folder, you will find the output of phenix.refine, and a folder called hydrate_structure. Inside this folder, you will find lilly_refined_final.pdb, which is the final structure of interest post refinement.

## Troubleshooting

Depending on your installed version of Rosetta, the `-restore_pre_talaris_2013_behavior true` flag in the spades_args.txt file in input_files can cause a segmentation fault when using rosetta_scripts. It is unknown why this occurs, but to be safe we have included the flag within the HydrateMover so you can remove it - although ideally you should keep it if you can for the density functions. This flag should still work regardless with the hydrate commandline function. 

## Contact

If you have any issues with the installation, demo or otherwise, please either submit a ticket above or contact me at lucas.rudden@epfl.ch



