#!/usr/bin/env bash

cp ATLAS_2017_I1512776_abhack-vn4.cc ATLAS_2017_I1512776.cc
if rivet-buildplugin ATLAS_2017_I1512776.cc -Wno-deprecated-declarations; then
echo "compiled!"
if [ $# -eq 0 ]
  then
    echo "No arguments supplied: running athena RivetAnalysis.py"
    athena RivetAnalysis.py
else
    echo "argument supplied: running athena RivetAnalysis_test1file.py"
    athena RivetAnalysis_test1file.py
fi

fi
