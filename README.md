# Stanford_Raytracer

Welcome to the Stanford Raytracer. This ray tracing program computes VLF ray paths in the magnetosphere.
This ray tracer has been moved to run on CU's supercomputer (Blanca via Summit -- see more info here).

I would strongly recommend following the instructions below to running this repo on the supercomputer. 
If you need to make changes to the source code (all Fortran based) then follow the instructions at the very end to
install this repo on your local machine and modify it. 

Firstly, this requires gfortran-8. A lot of new systems are using gfortran-9, so you might need to install the earlier version and make sure to update your priorities. Follow this link (change everywhere you see gcc to gfortran)
https://wiki.usask.ca/display/MESH/How-to+install+and+run+multiple+versions+of+GNU+compilers+in+an+Ubuntu+environment

You can check which version you are using by running in terminal: 
``` terminal
gfortran --version
```

First, make sure you have the repo. I'd reccommend forking it from this page into your own github account. 
Once you do that, in your own repo, find the url and enter into terminal (in whatever directory you want the project to build)
```terminal
git clone https://github.com/your-username/Stanford_Raytracer.git
``` 

Next, navigate inside the repo: 
``` terminal
cd Stanford_Raytracer
```

Next, run the following line to update the damping submodule:
``` terminal
git submodule init  
git submodule update  
``` 

The damping submodule also has submodules..
run:
``` terminal
cd damping
git submodule init  
git submodule update  
``` 

Finally, build:
``` terminal
cd ..
make
```

That's it! Reach out if something breaks...