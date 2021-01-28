# Stanford Raytracer

## Requirements

You must be using an x86 system to build this repo. It will not work on other hardware architectures.
Generally, if you're on a desktop or laptop computer, you're on an x86 system.

Mac/Windows:

* [Install Docker Desktop](https://www.docker.com/products/docker-desktop). Docker must be installed and *running*.

Linux:

* Docker Engine must be installed. Links to installation instructions for several distributions are given below:
  * [Ubuntu](https://docs.docker.com/engine/install/ubuntu/)
  * [Debian](https://docs.docker.com/engine/install/debian/)
  * [Fedora](https://docs.docker.com/engine/install/fedora/)

## Using This Repo

### Building and Running the Stanford Raytracer

If you just need to build and run the Standford Raytracer (or call it from another program, like a python script), do the following:

1. Ensure Docker is running  
linux - ``` docker info ```  
2. In a terminal, navigate to this project's directory and run:

    ```terminal
    docker build -t srt .
    ```

3. To run the Stanford Raytracer, run the command:

    ```terminal
    docker run srt {your arguments here}
    ```

    ```terminal
    docker run srt --outputper=1 --dt0=0.001 --dtmax=0.1 --tmax=1.5 --root=2 --fixedstep=0 --maxerr=0.0005 --modelnum=1 --maxsteps=2000 --minalt=6846200 --inputraysfile="/home/rileyannereid/workspace/SR_output/data_analysis_spring2021/2020-09-14 22:55:00/ray_inpfile.txt" --outputfile="/home/rileyannereid/workspace/SR_output/data_analysis_spring2021/2020-09-14 22:55:00/_ray_out_mode1.ray" --yearday=2020258 --milliseconds_day=82500000 --use_tsyganenko=1 --use_igrf=1 --tsyganenko_Pdyn=4 --tsyganenko_Dst=1 --tsyganenko_ByIMF=0 --tsyganenko_BzIMF=-5 --tsyganenko_W1=0.132 --tsyganenko_W2=0.303 --tsyganenko_W3=0.083 --tsyganenko_W4=0.07 --tsyganenko_W5=0.211 --tsyganenko_W6=0.308 --ngo_configfile="ngoconfig.in"
    ```
