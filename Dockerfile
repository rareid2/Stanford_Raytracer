FROM debian:buster AS RT-BUILDER

# this is where we will run our commands in this file
WORKDIR /usr/src/app

# Update packages and install and say yes to everything
RUN apt update && apt install -y \
    build-essential \ 
    gfortran \
    git

# Copy project over to the docker image
COPY Makefile .
COPY xform_double xform_double
COPY xform xform
COPY tsyganenko tsyganenko
COPY tricubic-for tricubic-for
COPY lapack-3.2.1 lapack-3.2.1
COPY iri2007 iri2007
COPY gcpm gcpm
COPY fortran fortran
RUN mkdir -p bin

# build the raytracer
RUN make

FROM python:3.7 AS XFORM-BUILDER

RUN apt update && apt install -y \
    gfortran

WORKDIR /usr/src/app

# move to submodule and build it
COPY python_setup/libxformd libxformd
WORKDIR /usr/src/app/libxformd/xform_double
RUN make && make shared

FROM python:3.7

RUN apt update && apt install -y \
    libgfortran5

WORKDIR /usr/src/app

# get the requirements txt for python
COPY python_setup/requirements.txt python_setup/requirements.txt

# get the needed python packages
RUN pip3 install -r python_setup/requirements.txt

COPY python_setup/*.py python_setup/

COPY --from=XFORM-BUILDER /usr/src/app/libxformd/libxformd.a python_setup/libxformd/libxformd.a
COPY --from=XFORM-BUILDER /usr/src/app/libxformd/libxformd.so python_setup/libxformd/libxformd.so
COPY --from=XFORM-BUILDER /usr/src/app/libxformd/xflib.py python_setup/libxformd/xflib.py

COPY --from=RT-BUILDER /usr/src/app/bin/getBfield bin/getBfield
COPY --from=RT-BUILDER /usr/src/app/bin/raytracer bin/raytracer

#CMD ["/bin/echo", "hI I aM raYTracEr"]
ENTRYPOINT ["python3", "/usr/src/app/python_setup/ray_input.py"]