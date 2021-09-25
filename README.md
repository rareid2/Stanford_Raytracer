# Stanford_Raytracer
Welcome to the Stanford Raytracer. This ray tracing program computes VLF ray paths in the magnetosphere.
This ray tracer has been moved to run on CU's supercomputer (Blanca and Summit -- see more info [here](https://curc.readthedocs.io/en/latest/index.html)).

I would strongly recommend following the instructions below to running this on the supercomputer. 
If you need to make changes to the source code (Fortran) then follow the instructions at the very end to
install this repo on your local machine and modify it. 

## First a quick tour of the folders in here 
- bin - executables (the docker image pulls the ./raytracer from here)
- fortran - all fortran source code (raytracer, plasmasphere model adapters, etc.)
- gcpm - GCPM (a plasmasphere model, the simplification to create mode 6 occurs in fortran source code)
- iri2007 - IRI (an ionosphere model, required by GCPM)
- lapack-3.2.1 - Lapack and BLAS, used for solving Ax=b
- matlab - a raytracer written in matlab. I think it works, but I've never used it...
- tricubic-for - Tricubic interpolator
- tsyganenko - Tsyganenko magnetosphere model
- xform - coordinate transformation tools (single precision)
- xform_double - coordinate transformation tools (double precision), this is also used inside of python_setup for coordinate transforms with a python wrapper... might be redudant but oh well

and the files:

- Makefile - this runs the make command to build the raytracer
- Dockerfile - buils the raytracer into an image
- .dockerignore, .gitignore, .gitmodules - docker and git specific files to keep it nice
- default_config.json - an example of a raytracer input file
- default_config_description.md - a description of every single parameter in more detail than you likely want for the config.json file
- exmaple_blanca_run.sh - an example job script to run on blanca nodes
- exmaple_summit_run.sh - an example job script to run on summit nodes

## Getting your rays ready to run

If you don't already have a supercomputing account, skip to step 3 to start that process, it can take a while. Then come back to step 1 while you wait for your account to be activated. 

1. *Create config file:* Let's start with the config file. It's easier to create this on your local machine and then copy it over. 

  - First, wherever you like to work on code, make a file called config.json.

  - This will have all the settings for your raytracer run. There is a [default_config.json](https://github.com/rareid2/Stanford_Raytracer/blob/master/default_config.json) in this repo that you can copy from and re-work as you like. You can copy over the default settings in the top of the file, and then edit the rays section. Details on these settings can be found here: [default_config_description.md](https://github.com/rareid2/Stanford_Raytracer/blob/master/default_config_description.md)

  - If you're running lots of rays, I'd reccommend generating this section with your language of choice. To generate the rays section from python, see [here](https://www.kite.com/python/answers/how-to-dump-a-dictionary-to-a-json-file-in-python). For matlab, see [here](https://www.mathworks.com/help/matlab/ref/jsonencode.html). I chose json format because it's compatible with a lot of languages and easy to edit from a text editor (and bc a software developer told me too).

  - Once you are happy with the config file, move on to the job script file. 

2. *Create job script:* If you're running on blanca, take a look at [example_blanca_run.sh](https://github.com/rareid2/Stanford_Raytracer/blob/master/example_blanca_run.sh) and for summit take a look at [example_summit_run.sh](https://github.com/rareid2/Stanford_Raytracer/blob/master/example_summit_run.sh).

  - Make any edits you see fit (name, date, title of run). Finally, make sure your config file is titled 'config.json'. If it's not, you're going to need to change the run command to the correct name of the config file and path. This run script assumes you're going to put it in the projects directory, so let's do that next. 

## Getting the raytracer on summit
3. *Get an account:* You need an account. Follow the instructions [here](https://rcamp.rc.colorado.edu/accounts/account-request/create/organization) to get one.

4. *File transfer:* Next, move the run script and config file to supercomputer. 

  - Run this command from terminal to move it (you'll need log in info) `scp <path_to_config> <path_to_job_script> <username>@login.rc.colorado.edu:/projects/<username>/` 

5. *Log in:* Now, actually log in to the supercomputer:

  - `ssh <username>@login.rc.colorado.edu`

  - Learn a little about the file structure [here](https://curc.readthedocs.io/en/latest/compute/filesystems.html). You have a home directory, projects directory, and a scratch directory. I'd reccommend installing the raytracer image into the projects directory and keeping your job scripts and config files in here too.

  - Make sure the files are where they are supposed to be: 

  - `cd /projects/<username>/`

  - `ls` 

  - If you see a config file here and a job script you're good to go!

6. *Load needed modules:* We now need to pull the latest image of the raytracer into summit. We have to do this from a compile node, so go ahead and run: 

  - `ssh scompile`

   - Note that when you run this, it will take you right back to your home directory. Let's move over to projects again:

  - `cd /projects/<username>/` 

  - Now let's get singularity module loaded (a containerization software supported by summit - more information [here](https://curc.readthedocs.io/en/latest/software/ContainerizationonSummit.html?highlight=singularity#singularity)):

  - `module load singularity/3.6.4`

7. *Pull latest image:* And finally, let's get the image (this might take a sec, it's about 1GB):

  - `singularity pull docker://rareid/raytracer:latest`

  - If you've got a lot of images, or you re-ran this command a lot, you may encounter an error that says you don't have enough space for this image. If that happens, try to purge your singularity cache:

  - `singularity cache clean` 

8. *Submit job:* That's everything we need! Let's submit the job. 

  - Type `cntrl+d` to exit out of the compile node. 

  - If you're running on blanca, you need to load that module real quick: 

  - `module load slurm/blanca` 

  - Now, run `sbatch <jobfile.sh>` to submit your job.

  - And you're done! Terminal output (ray exit code, etc.) will be in the `REPORT` file when the job finishes, in your project directory. The actual ray output will be in a folder created where you specificed in the config file (likely in your project directory as well). For some fun plotting tools, and output file parsing, check out this repo [here](https://github.com/rareid2/SR_interface).

## Running the raytracer locally using docker

*tbd*

## Building the raytracer locally -- do not reccommend

*in progress* 

1. Make sure you have gfortran-8 installed and update your priority list (decides which gfortran install to use). Here's some help on that here (change gcc to gfortran). 

2. Double check your default version is now gfortran-8: `gfortran --version` 

3. Clone this repo where you'd like it (`<path_to_raytracer>`). I reccommend a folder on your computer that you code in. If you have a github account or want to use github, you could instead fork this repo and clone from your own remote repo rather than this one. 

4. Next, head on in to the raytracer wherever you put it. `cd <path_to_raytracer>/Stanford_Raytracer` 

5. Next, we just need to build it. run `make clean` and `make` 