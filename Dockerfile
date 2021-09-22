FROM debian:buster

# this is where we will run our commands in this file
WORKDIR /usr/src/app

# Update packages and say yes to everything
RUN apt update -y

# Install GCC and say yes to everything
RUN apt install build-essential -y

# Verify that GCC is installed version 8
RUN gcc --version

# install gfortran
RUN apt install gfortran -y

# Snag git so we can install git-based pip packages (libxformd)
RUN apt install git -y

# now we need python3
RUN apt-get install python3 -y

# we need pip also
RUN apt install python3-pip -y

# Copy project over to the docker image
COPY . .

# get the needed python packages
RUN pip3 install -r python_setup/requirements.txt

# now we need to add in the git submodule for coordinate transforms
RUN git submodule init

# build!
RUN make

# default command? 
CMD ["echo", "Hello im the stanford raytracer"]