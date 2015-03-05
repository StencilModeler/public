#!/bin/bash

echo "========="
echo "== EXP =="
echo "========="
./model.1thr.1xCore.EXP_INTERP 64 64 64 16 16 16 20 1
echo
echo "========="
echo "== LIN =="
echo "========="
./model.1thr.1xCore.LIN_INTERP 64 64 64 16 16 16 20 1
echo
echo "========="
echo "== LOG =="
echo "========="
./model.1thr.1xCore.LOG_INTERP 64 64 64 16 16 16 20 1
echo
echo "========"
echo "== NO =="
echo "========"
./model.1thr.1xCore.NO_INTERP 64 64 64 16 16 16 20 1
echo
echo "=========="
echo "== NONO =="
echo "=========="
./model.1thr.1xCore.NO_INTERP.NO_SCALE 64 64 64 16 16 16 20 1
