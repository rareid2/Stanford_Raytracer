FROM debian:buster

WORKDIR /usr/src/app

# Update packages and say yes to everything
RUN apt update -y

# Install GCC and say yes to everything
RUN apt install build-essential -y

# Verify that GCC is installed
RUN gcc --version

# install gfortran
RUN apt install gfortran -y

# Snag git so we can install git-based pip packages (i.e. Spacepy fork)
RUN apt install git -y

# Copy project over to the docker image
COPY . .

# build!
RUN make

CMD ["/bin/bash"]