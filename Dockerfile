FROM debian:buster AS build

WORKDIR /usr/src/app

# a fix for apt-get to ensure it can properly run update
#RUN printf "deb http://archive.debian.org/debian/ buster main\ndeb-src http://archive.debian.org/buster/ buster main\ndeb http://security.debian.org buster/updates main\ndeb-src http://security.debian.org buster/updates main" > /etc/apt/sources.list

# Update packages and say yes to everything
RUN apt-get update -y && apt-get install -y \
    build-essential \
    gfortran

# make a bin folder for the raytracer binary
RUN mkdir -p bin

# copy everything from our local directory into the container's current directory
# the container's current directory was set at the WORKDIR instruction
COPY . .

# build!
RUN make

FROM debian:buster

RUN apt-get update && apt-get install -y \
    libgfortran4

WORKDIR /usr/src/app

COPY --from=build /usr/src/app/bin/raytracer .

ENTRYPOINT [ "/usr/src/app/raytracer" ]
