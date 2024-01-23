#!/bin/bash
sudo docker run -v ./container_data:/projects/$USER <dockerusername>/<tag> /projects/$USER/config.json