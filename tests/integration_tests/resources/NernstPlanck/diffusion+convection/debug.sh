#!/bin/bash

source /home/jyli/.bash_copss_midway2
source $LIBMESH_DIR/examples/run_common.sh

echo ------------------------------------------------------------------
gdb -ex=r --args $COPSS_DIR/src/copss-POINTPARTICLE-opt                  \
--disable-perflog                          \
-ksp_type preonly                           \
-pc_type lu                                 \
-pc_factor_mat_solver_package superlu_dist  \
-ksp_monitor -eps_monitor \
-memory_view -log_view 2>&1 | tee log_history.txt

