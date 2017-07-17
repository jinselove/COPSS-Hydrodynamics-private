#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh



# 1. *****************************************************************************************
# ========================== FieldSplit Schur Complement Solver ====================
# with schur_user_ksp=true
# mesh          :  50x10x50                                 100x20x100(5+M dofs)
# solver time   :  82.4527 s                                4520.9760 s
# iteration step:  42 KSP Residual norm 1.343672186970e-06  44 KSP Residual norm 2.435167699399e-06
# This works very fast only for schur_user_ksp = true, but needs more iterative steps to converge;
# if false, the iteration is less for convergence, but total time used is more!
# For some special mesh configuration, e.g. 40X16X16 or coarser, it fails to converge! (if true)
#-malloc_dump # outputs all unfreed memory after PetscFinalize() to see the leak.
#-memory_view -log_view

# echo ------------------------------------------------------------------
# echo ------------------------ Use 4 processors ------------------------
# mpirun -n 1 --bind-to core ./example-opt \
# -ksp_type gmres  -ksp_gmres_restart 100   \
# -pc_type fieldsplit                       \
# -pc_fieldsplit_type schur                 \
# -pc_fieldsplit_schur_fact_type lower      \
# -fieldsplit_0_ksp_type cg            \
# -fieldsplit_0_pc_type bjacobi               \
# -fieldsplit_1_pc_type bjacobi             \
# -fieldsplit_1_ksp_type cg              \
# -malloc_dump -malloc_debug -ksp_monitor -fp_trap \
# -memory_view -log_view 2>&1 | tee log_history.txt
###-start_in_debugger noxterm
###-fieldsplit_0_ksp_monitor
###-fieldsplit_0_ksp_max_it 100
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./example-opt 


# TEST: **************************************************************************************
# use lu exact solve for velocity(A00*u=f), and reuse the LU matrices
# NOTE: for fieldsplit_1, we solve Mp*p = g, which is not exact!!!
#echo ------------------------------------------------------------------
#echo ------------------------ Use 4 processors ------------------------
#mpiexec -n 4 -bind-to hwthread -map-by socket:1 ./example-opt  \
#-ksp_type gmres  -ksp_gmres_restart 100   \
#-pc_type fieldsplit                       \
#-pc_fieldsplit_type schur                 \
#-pc_fieldsplit_schur_fact_type full      \
#-fieldsplit_0_ksp_type preonly            \
#-fieldsplit_0_pc_type lu                  \
#-fieldsplit_0_pc_factor_mat_solver_package superlu_dist \
#-fieldsplit_1_ksp_type preonly            \
#-fieldsplit_1_pc_type lu                  \
#-fieldsplit_1_pc_factor_mat_solver_package superlu_dist \
#-malloc_dump -malloc_debug -ksp_monitor \
#-memory_view -log_view 2>&1 | tee log_history.txt


# 2. *****************************************************************************************
# Direct solver with super_LU_dist, note this is only for coarse mesh!
# When this solver is used, go to xxx_control.in file and turn
# user_defined_pc = false, and schur_complement = false
#-mat_superlu_dist_colperm MMD_AT_PLUS_A     \ Increase computational time
#-mat_superlu_dist_statprint - print factorization information

echo ------------------------------------------------------------------
echo ------------------------ Use 1 processors ------------------------
mpirun -n 50 ./example-opt                   \
--disable-perflog                          \
-ksp_type preonly                           \
-pc_type lu                                 \
-pc_factor_mat_solver_package superlu_dist  \
-ksp_monitor -eps_monitor \
-memory_view -log_view 2>&1 | tee log_history.txt
###-info -log_summary -ksp_view_pmat binary -fp_trap
### $PETSC_DIR/$PETSC_ARCH/bin/


### lldb is better than gdb un MacOS
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./example-opt \
#-ksp_type preonly                           \
#-pc_type lu                                 \
#-pc_factor_mat_solver_package mumps  \
#-ksp_monitor 2>&1 | tee log_history.txt
###-log_summary
#run -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package superlu_dist


#  -------------------------------------------------------------------------------------------
# the followings are for test of solvers.
#  -------------------------------------------------------------------------------------------
# this is Super SLOW! However, we tested 20x8x8 case, it converges perfectly(schur_user_ksp = false)
#mpirun -np 6 ./example-opt -ksp_type gmres -ksp_gmres_restart 100 -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type bjacobi -fieldsplit_1_pc_type bjacobi -fieldsplit_1_ksp_type gmres -fieldsplit_1_inner_ksp_type preonly -fieldsplit_1_inner_pc_type bjacobi -ksp_monitor
#0 KSP Residual norm 5.577283041725e-01
#1 KSP Residual norm 3.454666793318e-01
#2 KSP Residual norm 2.686535098544e-01
#3 KSP Residual norm 1.699660239329e-01
#4 KSP Residual norm 1.128272384313e-01
#5 KSP Residual norm 7.284399151258e-02
#6 KSP Residual norm 4.195003239102e-02
#7 KSP Residual norm 2.081072032792e-02
#8 KSP Residual norm 1.190515122415e-02
#9 KSP Residual norm 8.271962626272e-03
#10 KSP Residual norm 4.124557343271e-03
#11 KSP Residual norm 2.333149219441e-03
#12 KSP Residual norm 1.433337885639e-03
#13 KSP Residual norm 7.377561711995e-04
#14 KSP Residual norm 4.207075229197e-04
#15 KSP Residual norm 1.980294252024e-04
#16 KSP Residual norm 9.292470721883e-05
#17 KSP Residual norm 4.275893946337e-05
#18 KSP Residual norm 2.074274111312e-05
#19 KSP Residual norm 9.108299154960e-06
#20 KSP Residual norm 4.650248554765e-06
#21 KSP Residual norm 2.047444308323e-06
#22 KSP Residual norm 9.115176809842e-07
#23 KSP Residual norm 4.306531645198e-07


# From SNES example 70, Super SLOW!
#mpiexec -np 4 ./example-opt -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type bjacobi -fieldsplit_1_pc_type jacobi -fieldsplit_1_inner_ksp_type preonly -fieldsplit_1_inner_pc_type jacobi -fieldsplit_1_upper_ksp_type preonly -fieldsplit_1_upper_pc_type jacobi



# ------------------------------------- test of PCFieldSplit -----------------------------------
# test: 1.0
# Cohouet & Chabard, Some fast 3D finite element solvers for the generalized Stokes problem, 1988.
# Modification: ml -> hypre
# This works for schur_user_ksp = false or true,
# but slower than pc_fieldsplit_type = schur/multiplicative
#mpirun -n 4 ./example-opt    \
#-ksp_type fgmres -ksp_gmres_restart 100  \
#-pc_type fieldsplit           \
#-pc_fieldsplit_type multiplicative  \
#-fieldsplit_0_pc_type hypre      \
#-fieldsplit_0_ksp_type preonly  \
#-fieldsplit_1_pc_type jacobi    \
#-fieldsplit_1_ksp_type preonly  \
#-fieldsplit_1_ksp_monitor_true_residual \
#-fieldsplit_1_ksp_converged_reason \
#-ksp_monitor_true_residual #-ksp_view -ksp_converged_reason


# test: 2.0
# Elman, Multigrid and Krylov subspace methods for the discrete Stokes equations, 1994.
# If schur_user_ksp = false, it converges fast, but the solution in incorrect.
#mpirun -n 4 ./example-opt    \
#-ksp_type gmres -ksp_gmres_restart 100  \
#-pc_type fieldsplit                     \
#-pc_fieldsplit_type multiplicative      \
#-fieldsplit_0_pc_type hypre             \
#-fieldsplit_0_ksp_type preonly          \
#-fieldsplit_1_pc_type jacobi            \
#-fieldsplit_1_ksp_type cg          \
#-ksp_monitor #-ksp_view -ksp_converged_reason


# -------------------------------- Comments: --------------------------------
# 2.0 is basically a variant of 1.0 by simply replace additive by multiplicative.
# Similarly, we can change other options:
# -pc_fieldsplit_type additive/multiplicative/schur
# -fieldsplit_0_pc_type hypre/ml/gamg (if use gamg, -fieldsplit_0_pc_gamg_sym_graph true)
# ---------------------------------------------------------------------------


# test: 3.1
# May and Moresi, Preconditioned iterative methods for Stokes flow problems arising in
# computational geodynamics, 2008.
# Olshanskii, Peters, and Reusken, Uniform preconditioners for a parameter dependent saddle point
# problem with application to generalized Stokes interface equations, 2006.
# Modification: ml -> hypre; minres -> gmres. This is slow comparing with others
#mpirun -np 4 ./example-opt      \
#-ksp_type gmres -ksp_gmres_restart 100  \
#-pc_type fieldsplit             \
#-pc_fieldsplit_type schur       \
#-fieldsplit_0_pc_type hypre      \
#-fieldsplit_0_ksp_type preonly  \
#-fieldsplit_1_user_pc none      \
#-fieldsplit_1_ksp_type gmres   \
#-pc_fieldsplit_schur_fact_type diag  \
#-ksp_monitor #-ksp_view -ksp_converged_reason


# test: 3.2
# May and Moresi, Preconditioned iterative methods for Stokes flow problems arising in
# computational geodynamics, 2008.
# Modification: -fieldsplit_0_pc_gamg_sym_graph true; minres->gmres; 1_pc_type none->jacobi (slow)
#mpirun -np 4 ./example-opt      \
#-ksp_type gmres -ksp_gmres_restart 100  \
#-pc_type fieldsplit             \
#-pc_fieldsplit_type schur       \
#-fieldsplit_0_pc_type gamg      \
#-fieldsplit_0_pc_gamg_sym_graph true   \
#-fieldsplit_0_ksp_type preonly  \
#-fieldsplit_1_pc_type jacobi      \
#-fieldsplit_1_ksp_type gmres   \
#-pc_fieldsplit_schur_fact_type lower  \
#-ksp_monitor #-ksp_view -ksp_converged_reason


# test: 3.3
# May and Moresi, Preconditioned iterative methods for Stokes flow problems arising in
# computational geodynamics, 2008.
#mpirun -np 4 ./example-opt -ksp_type gmres -ksp_gmres_restart 100 -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_pc_type gamg -fieldsplit_0_ksp_type preonly -fieldsplit_1_pc_type none -fieldsplit_1_ksp_type minres -pc_fieldsplit_schur_factorization_type upper -pc_fieldsplit_detect_saddle_point -ksp_view -ksp_monitor -ksp_converged_reason


# test: 3.4
# May and Moresi, Preconditioned iterative methods for Stokes flow problems arising in
# computational geodynamics, 2008.
# Kay, Loghin and Wathen, A Preconditioner for the Steady-State N-S Equations, 2002.
# Elman, Howle, Shadid, Shuttleworth, and Tuminaro, Block preconditioners based on approximate
# commutators, 2006.
# The same as test 3.1, this is also extremely slow due to fieldsplit_1_pc_type lsc !
#mpirun -np 4 ./example-opt -ksp_type gmres -ksp_gmres_restart 100 -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_pc_type gamg -fieldsplit_0_pc_gamg_sym_graph true -fieldsplit_0_ksp_type preonly -fieldsplit_1_pc_type lsc -fieldsplit_1_ksp_type minres -pc_fieldsplit_schur_factorization_type upper -pc_fieldsplit_detect_saddle_point -ksp_view -ksp_monitor -ksp_converged_reason

