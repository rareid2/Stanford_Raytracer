#!/bin/bash
#SBATCH --nodes=1                               # Number of requested nodes (keep this at 1 for raytracer bc i couldnt figure out cross node computing)
#SBATCH --time=23:59:00                         # Max wall time (maybe increase this for over 500,000 rays)
#SBATCH --qos=blanca-lair                       # Specify testing QOS (this is for lair blanca nodes)
#SBATCH --ntasks=40                             # Number of tasks per job (leave this at 40, the raytracer will use them all)
#SBATCH --job-name=raytracer                    # Job submission name (whatever you like!)
#SBATCH --output=REPORT.%j.out                  # Output file name with Job ID (this is the ray output, you should save and look at this, it shows the exit codes of each ray)

# Written by:   Riley Reid
# Date:         22 Sept 2021
# Purpose:      run the stanford raytracer on blanca

# load singularity (needed to run this image)   
module load singularity/3.6.4

# go to projects directory! where this image and the config.json should be
cd /projects/$USER/
# execute the image using the python script which reads the config file and calls the raytracer binary 
singularity exec raytracer_latest.sif config.json