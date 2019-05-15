export RIVET_ANALYSIS_PATH=$PWD
rivet-mkhtml Acer.yoda
rm ~/temp_dump_for_cpin/*
cp rivet-plots/ATLAS_2017_I1512776/*.png ~/temp_dump_for_cpin/
cp rivet-plots/*.html ~/temp_dump_for_cpin/
