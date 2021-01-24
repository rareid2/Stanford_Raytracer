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
2. In a terminal, navigate to this project's directory and run:

    ```terminal
    docker build -t srt .
    ```

3. To run the Stanford Raytracer, run the command:

    ```terminal
    docker run srt {your arguments here}
    ```
