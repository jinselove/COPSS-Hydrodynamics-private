.. _getting-started:

Getting Started
================

Installation
--------------
**COPSS** is written on [LIBMESH](http://libmesh.github.io/) framework. It also requires [PETSc](https://www.mcs.anl.gov/petsc/index.html) for parallel linear equation solvers and [SLEPc](http://slepc.upv.es/) for scalable Eigenvalue computations. Before installing **COPSS**, you need to install **LIBMESH** together with **PETSc** and **SLEPc**. To achieve the best parallel performance of COPSS, we suggest install it on a Linux cluster environment.

0. System environment prep

 - [CMAKE](https://cmake.org/) (e.g., but not necessarily, version 3.6.2)
 - [GCC](https://gcc.gnu.org/) (e.g., but not necessarily, version 6.2)
 - [PYTHON](https://www.python.org/) (python 2)
 - [OPENMPI](https://www.open-mpi.org/) (e.g., but not necessarily, version 2.0.1)


1. Install PETSc

 - Download PETSC's latest release ( version 3.7.4 or later ) from [PETSc download](https://www.mcs.anl.gov/petsc/download/index.html), or git clone PETSc repository:
	 
        `mkdir $HOME/projects/`

        `cd $HOME/projects/`

        `git clone -b maint https://bitbucket.org/petsc/petsc petsc`
     
 - Configure PETSc:
 
        `cd $HOME/projects/petsc`

        `./configure --with-cc=mpicc --with-cxx=mpicxx --with-mpiexec=mpiexec --with-fc=mpif90 �download-fblaslapack --download-scalapack --download-mumps --download-superlu_dist --download-hypre --download-ml --download-parmetis --download-metis --download-triangle --download-chaco --with-debugging=0`

        **And then follow the instructions on screen to install and test the package.**
        
 - Export environment variables:
	
        `export PETSC_DIR=*/path/to/PETSC*`

        `export PETSC_ARCH=*PETSC_ARCH_NAME*`
	
        **Add the above codes to ~/.bashrc and `source ~/.bashrc` before next step. (`*/path/to/PETSC*` and `*PETSC_ARCH_NAME*` can be found on the screen after installation.)**

 - If you meet any trouble, please refer to [PETSC installation](https://www.mcs.anl.gov/petsc/documentation/installation.html).
	
2. Install SLEPc
 - Download SLEPC's latest release (version 3.7.3 or later) from [SLEPc download](http://slepc.upv.es/download/download.htm), or git clone PETSc repository:
	 
        `cd $HOME/projects/`


        `git clone -b maint https://bitbucket.org/slepc/slepc slepc`
	 
 - Configure PETSc:
 
        `cd $HOME/projects/slepc`

        `./configure`

        **And then follow the instructions on screen to install the package**
	  
 - Export environment variables:
	
        `export SLEPC_DIR=*/path/to/SLEPC*`

        `export SLEPC_ARCH=*SLEPC_ARCH_NAME*`

        **Add the above codes to ~/.bashrc and `source ~/.bashrc` before next step. (`*/path/to/SLEPC*` and `*SLEPC_ARCH_NAME*` can be found on the screen after installation.)**

 - Test the package (not necessary but recommended)
 
        `make test`

If you meet any trouble, please refer to [SLEPC installation](http://slepc.upv.es/documentation/instal.htm).

3. Install LIBMESH

 - Download LIBMESH's latest release ( version 1.1.0 or later ) from [LIBMESH download](https://github.com/libMesh/libmesh/releases), or git clone PETSc repository:
	 
        `cd $HOME/projects/`

        `git clone git://github.com/libMesh/libmesh.git`
	 
 - Build LIBMESH:
 
        `cd $HOME/projects/libmesh`

        `./configure -prefix=$HOME/projects/libmesh/libmesh-opt --enable-optional --enable-vtk  --enable-gzstream --enable-trilinos --disable-strict-lgpl --enable-laspack --enable-capnproto --enable-trilinos --enable-nodeconstraint --enable-perflog --enable-ifem --enable-petsc --enable-blocked-storage --enable-slepc --enable-unique-id --enable-unique-ptr --enable-parmesh 2>&1  | tee my_config_output_opt.txt`

        (Read the configuration output, make sure **PETSC** and **SLEPC** is enabled).

        **And then follow the instructions on screen to install and test the package.**

 - Export environment variables:
	
        `export LIBMESH_DIR=*/path/to/LIBMESH*`

        **Add the above codes to ~/.bashrc and `source ~/.bashrc` before next step. (`*/path/to/PETSC*`can be found on the screen after installation.)**

 - If you meet any trouble, please refer to [LIBMESH installation](https://libmesh.github.io/installation.html), or reach out to **LIBMESH** community for help.

4. Install COPSS-hydrodynamics

 - Download the latest COPSS codes
 
        `cd /path/to/where/you/want/to/install/copss`

        `git clone https://bitbucket.org/COPSS/copss-hydrodynamics-public.git`

 * **On Master branch (only point particle systems are available)**


    - Compile the codes (on Master branch)

        `cd /path/to/copss/src/`

        `make`



    - Run an diffusion example of beads

        `cd /path/to/copss/examples/general_point_particle/beads/`

        `cp /path/to/copss/src/example-opt $PWD`

        `cp /path/to/copss/tools/run.sh $PWD`

        `bash run.sh` (You can define how many cores you want to run on in **run.sh**)

	
	You need to set up your system in **point_particle_control.in** and **datafile (e.g., point_particle_data.in)**. More details can be found in our documentation.


 * **On Development branch (both PointParticle and RigidParticle are available)**
    - Compile the codes
 
        `cd /path/to/copss/src/`

        `make package=POINTPARICLE (for point particle systems)  Or make package=RIGIDPARTICLE (for rigid particle systems)`


    - Run an sedimentation example of rigid particles
 
        `cd /path/to/copss/examples/general_rigid_particle/sedimentation_benchmark/`

        `modify run.sh dependending on your system`

        `bash run.sh`
	  

**Build documentation**
-------------------------------------------
After you have build **COPSS-Hydrodynamics** you can further build the documentation
in `docs/` directory.

1. Doxygen
    The documentation built using **Doxygen** gives the code-leve details, including
    the code structure, class inheritance, details of functions, etc. To compile the documentation,
    make sure you have [Doxygen](http://www.stack.nl/~dimitri/doxygen/) ready, then:
    
        `cd [path-to/copss]/docs/doxygen/`
    
        `doxygen Doxyfile.bak`

    then you can view the documentation in an IE browser:

        `google-chrome [path-to-copss]/docs/doxygen/html/index.html`

2. Sphinx
   The documentation built using **Sphinx** gives the tutorial-level details, including
   features of the package, how to run a simulation, how to use a tool, etc. To compile Sphinx,
   make sure you have [Sphinx] (http://www.sphinx-doc.org/en/master/) ready, then:
   
    `cd [path-to-copss]/docs/sphinx`
   
    `make html`
   
   then you can view the documentation in an IE browser:
   
    `google-chrome [path-to-copss]/docs/sphinx/build/html/index.html`
    
   To modify Sphinx documentation, one needs to write the documentation in one of the rst files,
   for example, if the documentation is about how to run a simulation, it should be written into:
   
    `[path-to-copss]/docs/sphinx/source/tutorials.rst`
