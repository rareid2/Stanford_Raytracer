FROM debian:buster AS build

WORKDIR /usr/src/app

# Update packages and say yes to everything
RUN apt update -y && apt install -y \
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
    libgfortran-8-dev

WORKDIR /usr/src/app

COPY --from=build /usr/src/app/bin/raytracer .

ENTRYPOINT [ "/usr/src/app/raytracer" ]
