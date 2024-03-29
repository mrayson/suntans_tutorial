# Instructions to run SUNTANS on Pawsey

Steps to go through:

- Logging in to Pawsey
- Loading modules
- Compiling code
- Submitting a job to the queue using slurm
- Downloading data

General instructions for running on Pawsey (Magnus) are here:

- https://support.pawsey.org.au/documentation/display/US/Supercomputing+Documentation
- https://support.pawsey.org.au/documentation/display/US/Compiling+on+Magnus+and+Galaxy
- Cheat sheet: https://support.pawsey.org.au/documentation/download/attachments/34017103/Magnus_Quick_Reference_Guide.pdf?version=1&modificationDate=1536029582214&api=v2 

## Step 1) Logging in

Option 1: MobaXterm

Option 2: Use `ssh` from a linux terminal

`ssh username@setonix.pawsey.org.au`

## Step 2) Download the source code and this tutorial

Naviagate to your software folder: `cd $MYSOFTWARE` 

Git clone the code (suntans code and tutorial)

`git clone https://github.com/ofringer/suntans.git`

`git checkout variable_coriolis`



## Step 3) Loading modules

Use commands:

Set to the gnu programming environment. This is default on setonix.

`module list` list all currently loaded modules

`module avail` : list all available modules

Load netcdf: `module load netcdf-c/4.8.1`



### Download and build parmetis:

This installation is manual but you will learn how software packages are installed.

- ~~Download Parmetis: `wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/OLD/ParMetis-3.2.0.tar.gz`~~
- Download Parametis (backup):  `curl https://cloudstor.aarnet.edu.au/plus/s/Y9eVZuaKclvsY6W/download -o ParMetis-3.2.0.tar.gz`
- Untar the folder: `tar -xvf ParMetis-3.2.0.tar.gz`
- Change directory to the folder: `cd ParMetis-3.2.0`
- Change `mpicc` to `cc` in <parmetis-folder>/Makefile.in
- Compile the library: `make`
- done.

## Compiling suntans

- Download theexamples: `git clone https://github.com/mrayson/suntans_tutorial.git`
- Replace `suntans/main/Makefile` with `iwaves_shelf/Makefile.pawsey` i.e. `cp <path-to-suntans-tutorial>/iwaves_shelf/Makefile.pawsey <path-to-suntans-code/main>/Makefile`

- Update Makefile.in to:

```
#MPIHOME=
PARMETISHOME=/your/path/to/ParMetis-3.2.0
TRIANGLEHOME=
#NETCDF4HOME=
```

- Navigate to the `iwaves_shelf` test case folder
- Update `Makefile` $SUNTANSHOME variable e.g.
        `SUNTANSHOME='/home/mrayson/software/code/suntans/main`

- Run `make`

## Run an example (no python)

Clone the tutorial...

Go to the iwaves_shelf example `cd suntans_tutorial iwaves_shelf`

Edit the `Makefile`:

 - Replace `mpicc` with `cc`
 - Set the SUNTANSHOME variable to the path where SUNTANS is located

Update `iwaves_shelf_slurm.sh`

 - Set the SUNTANSHOME variable

Run the example:
  `sbatch iwaves_shelf_slurm.sh`
  
Check the job is running:
  `squeue -u mrayson`

Check the status of the job:
  `watch tail <name_of_log_file_in_slurm_script>`
  
The job will fail because your python environment is not setup correctly...

## Setting up python on Pawsey

Use a singularity (docker-like) container -- see the `iwaves_shelf/python_singularity` script












