#!/bin/bash

#set -x

rm .depend
rm *.o #*.m *.gmv #*.e
rm example-opt
rm -r example-opt.dSYM
rm example-dbg
rm -r example-dbg.dSYM

rm traceout* log_history.txt msd.txt
rm *velocity_profile*

rm output_test*
rm output*.e
rm output.e
rm pm_system*.e
rm particle_surface_mesh.*

rm *.vtu*
rm *.vtk
rm *.pvtu*
rm *.csv
rm *.txt
rm *.dat
rm out.*

rm binaryoutput*
rm vector_RIN.* vector_ROUT.*
rm *.e-s.*
