# Stanford_Raytracer

Instructions for building:  
First, clone the repo 
```terminal
git clone url
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