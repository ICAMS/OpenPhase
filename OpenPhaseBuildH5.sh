#!/bin/bash
git clone https://github.com/leethomason/tinyxml2.git
git clone https://github.com/BlueBrain/HighFive.git
git clone https://github.com/HDFGroup/hdf5.git --branch hdf5-1_12_1
cp tinyxml2/tinyxml2.h include/tinyxml2.h
cp tinyxml2/tinyxml2.cpp src/tinyxml2.cpp
cd hdf5
./configure
make
make install
cd ..
make clean
make SETTINGS="generic H5 shared" clean
make SETTINGS="generic H5 shared"
