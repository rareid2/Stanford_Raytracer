# Stanford_Raytracer
Welcome to the Stanford Raytracer. This ray tracing program computes VLF ray paths in the magnetosphere.
This ray tracer has been moved to run on CU's supercomputer (Blanca via Summit -- see more info here).

I would strongly recommend following the instructions below to running this on the supercomputer. 
If you need to make changes to the source code (Fortran) then follow the instructions at the very end to
install this repo on your local machine and modify it. 

## First a quick tour of the folders in here 
bin - executables (the docker image pulls the ./raytracer from here)

fortran - all fortran source code (raytracer, plasmasphere model adapters, etc.)

gcpm - GCPM (a plasmasphere model, the simplification to create mode 6 occurs in fortran source code)

iri2007 - IRI (an ionosphere model, required by GCPM)

lapack-3.2.1 - Lapack and BLAS, used for solving Ax=b

matlab - a raytracer written in matlab. i think it works, but I've never used it...

tricubic-for - Tricubic interpolator

tsyganenko - Tsyganenko magnetosphere model

xform - coordinate transformation tools (single precision)

xform_double - coordinate transformation tools (double precision), this is also used inside of python_setup for coordinate transforms with a python wrapper... might be redudant but oh well

and the files:

Makefile - this runs the make command to build the raytracer

Dockerfile - buils the raytracer into an image

.dockerignore, .gitignore, .gitmodules - docker and git specific files to keep it nice

config.json - an example of a raytracer input file

raytracer_settings.txt - a description of every single parameter in more detail than you likely want for the config.json file


## Getting the raytracer on summit
1. First, you need an account. Follow the instructions [here](https://rcamp.rc.colorado.edu/accounts/account-request/create/organization) to get one.
If you have on, log in from terminal, 
`ssh username@login.rc.colorado.edu`

Learn a little about the file structure [here](https://curc.readthedocs.io/en/latest/compute/filesystems.html). You have a home directory, projects directory, and a scratch directory. I'd reccommend installing the raytracer image into the projects directory, and keeping your job scripts and config files in here too. 

2. We need to pull the latest image of the raytracer into summit. We have to do this from a compile node, so go ahead and run: 
`ssh scompile`

Note that when you run this, it will take you right back to your home directory. Let's move over to projects:
`cd /projects/username` 

Now let's get singularity module loaded (a containerization software supported by summit - more information [here](https://curc.readthedocs.io/en/latest/software/ContainerizationonSummit.html?highlight=singularity#singularity)):
`module load singularity/3.6.4`

And finally, let's get the image (this might take a sec, it's about 1GB):
`singularity pull docker://rareid/raytracer:latest`

If you've got a lot of images, or you re-ran this command a lot, you may encounter an error that says you don't have enough space for this image. If that happens, try to purge your singularity cache:
`singularity cache clean` 

3. Now we have the raytracer image, and we are logged into summit. Our last two steps are to set up the config file and the job script.

Let's start with the config file. It's easier to create this on your local machine and then copy it over. 

First, wherever you like to work on code, make a file called `config.json.` 

This will have all the settings for your raytracer run. There is a `default_config.json` in this repo that you can copy from and re-work as you like. You can copy over the default settings in the top of the file, and then edit the rays section as you please. 

If you're running lots of rays, I'd reccommend generating this part with your coding language of choice. To generate the rays section from python, see [here](https://www.kite.com/python/answers/how-to-dump-a-dictionary-to-a-json-file-in-python). For matlab, see [here](https://www.mathworks.com/help/matlab/ref/jsonencode.html). I chose json because it's compatible with a lot of languages and easy to edit from a text editor (and bc a software developer told me too).

For details 


The job script sumbits to blanca or summit. There's an example for both blanca and summit included in this repo. 
To do that, first download the file `raytracer_job.sh`. I'd reccommend editing locally before moving to the supercomputer. Open any text editor and open the recently downloaded `raytracer_job.sh`.


For testing locally: (just for riley)
`docker build -t rareid/raytracer:latest` and `docker push rareid/raytracer:latest`


## Running the raytracer locally using docker

## Building the raytracer locally -- do not reccommend
1. Make sure you have gfortran-8 installed and update your priority list (decides which gfortran install to use). Here's some help on that here (change gcc to gfortran). 

2. Double check your default version is now gfortran-8: `gfortran --version` 

3. Clone this repo where you'd like it (`<path_to_raytracer>`). I reccommend a folder on your computer that you code in. If you have a github account or want to use github, you could instead fork this repo and clone from your own remote repo rather than this one. 

4. Next, head on in to the raytracer wherever you put it. `cd <path_to_raytracer>/Stanford_Raytracer` 

5. Next, we just need to build it. run `make clean` and `make` 