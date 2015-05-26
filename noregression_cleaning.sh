#!/bin/bash

make SANDY=1 THREADS=1 THRCORE=1 DEBUG=1

cd bin

./run.sh &> ./tmp.dat

diff ./baseline7.dat ./tmp.dat

echo
if [ $? -ne 0 ]; then echo "Non-regression tests: FAILED"
else echo "Non-regression tests: PASS"
fi

echo
echo "Done."
echo
