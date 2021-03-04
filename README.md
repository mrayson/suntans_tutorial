# Tutorial for setting up a SUNTANS test case

This assumes you have a linux-like environment. You could try on a Mac. Windows -- don't bother.

The basic steps are:
 
  - Install the relevant libraries, including python
  - Checkout SUNTANS source code from github (switch branches)
  - Compile SUNTANS
  - Run a test case
  - Inspect the model results visually

## Install C compiler and libraries

Need the following:
  
  - C compiler
  - MPI libraries
  - Make utility
  - NetCDF library
  - Parmetis

### Run the following in ubuntu to install:

`sudo apt install gcc make libnetcdf-dev libmpich-dev`

### Parmetis installation

This installation is manual but you will learn how software packages are installed.

- Download Parmetis: `wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/OLD/ParMetis-3.2.0.tar.gz`
- Untar the folder: `tar -xvf ParMetis-3.2.0.tar.gz`
- Change directory to the folder: `cd ParMetis-3.2.0`
- Compile the library: `make`
- done.

## Install Python and relevant libraries

- You will need `git` and `pip` installed first: `sudo apt install git python3-pip`
- Download and install a conda python installation. I like miniconda as its lightweight. `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
- Install the package by running the shell script e.g. `sh Miniconda<version>.sh`
- Create a new conda environment: `conda create --name sfoda`
- Activate the new environment: `conda activate sfoda`
- Install my `SFODA` python library that contains all of the relevant code for interfacing with suntans as well as major python libraries (numpy, scipy, xarray, etc): `pip install git+https://github.com/mrayson/sfoda.git`

## Checkout the SUNTANS source code

 - Navigate to the folder where you want to store the code: `cd suntans`
 - Clone the code from github: `git clone https://github.com/ofringer/suntans.git`
 - Switch to my branch: `git checkout variable_coriolis`

## Compile SUNTANS

 - Navigate to the `main` folder in the suntans code folder
 - Edit the file `Makefile.in` to point the compiler to the relevant library paths. For the installation described above, change it to:

```
MPIHOME=/usr
PARMETISHOME=/home/suntans/code/ParMetis-3.2.0
TRIANGLEHOME=
NETCDF4HOME=/usr
```
 
 - Compile the code by typing: `make`

This creates an executable file called `sun`. You can try running it by typing: `./sun`
 
It will return an error because you have not specfied the relevant inputs. This is next.

## Running a SUNTANS test case

- Download these examples: `git clone https://github.com/mrayson/suntans_tutorial.git`
- Navigate to the `tidal_front` example: `cd suntans_tutorial/tidal_front`
- Edit the Makefile. Point the first line to the `suntans/main` folder


---

Matt Rayson

University of Western Australia

February 2021






