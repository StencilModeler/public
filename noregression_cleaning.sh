#!/bin/bash

make SANDY=1 THREADS=1 THRCORE=1 DEBUG=1

cd bin

./run.sh &> ./tmp.dat

diff ./baseline5.dat ./tmp.dat
