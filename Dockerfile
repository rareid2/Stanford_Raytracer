FROM python:3.7-slim-buster

WORKDIR /usr/src/app

# Update packages and say yes to everything
RUN apt update -y

# Install GCC and say yes to everything
RUN apt install build-essential -y

# Verify that GCC is installed
RUN gcc --version

# install gfortran
RUN apt install gfortran -y

# Install pipenv so we can pull from the Pipfile
RUN pip install pipenv

# Snag git so we can install git-based pip packages (i.e. Spacepy fork)
RUN apt install git -y

# Copy project over to the docker image
COPY . .

# Download required python packages from the Pipfile
RUN pipenv update

# build!
RUN make

CMD ["sh", "-c", "python ./example_scripts/launch_multiplerays.py"]