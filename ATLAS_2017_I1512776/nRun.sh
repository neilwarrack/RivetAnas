export RIVET_ANALYSIS_PATH=$PWD

if [ "$1" == "-b" ]; then
    echo "attempting to build..."
    rivet-buildplugin ATLAS_2017_I1512776.cc -L. -lmm -lgsl{,cblas}
fi

if [ "$1" == "-t" ]; then
    rivet -a ATLAS_2017_I1512776 8Tev_noMpi_test_pythia8_events.hepmc
fi

if [ "$1" == "-h" ]; then
    rivet -a ATLAS_2017_I1512776 $2
fi
