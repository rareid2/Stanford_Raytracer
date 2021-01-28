# Stanford_Raytracer

Instructions for building:  
First, clone the repo 
```terminal
git clone url
``` 

Next, inside the repo, run the following to update the damping submodule

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