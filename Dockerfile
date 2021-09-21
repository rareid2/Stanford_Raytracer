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
RUN apt-get install python3

RUN pip install -r requirements.txt

# Copy project over to the docker image
COPY . .

# now we need to add in the git submodule
RUN git submodule init



# build!
#RUN make
RUN python
#CMD ["/bin/bash"]
