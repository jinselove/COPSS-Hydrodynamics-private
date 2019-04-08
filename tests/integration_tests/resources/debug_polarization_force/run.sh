#!/bin/bash

mpirun -n 28 $COPSS_DIR/src/copss-RIGIDPARTICLE-opt                  \
--disable-perflog                          \
-ksp_type preonly                           \
-pc_type lu                                 \
-pc_factor_mat_solver_package superlu_dist  \
-ksp_monitor -eps_monitor \
-memory_view -log_view 2>&1 | tee log_history.txt

#--disable-perflog                          \
