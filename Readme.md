**COPSS**
====================

**COPSS** (Continuum Particle Simulation Suite) is an open source,  [**LIBMESH**](http://libmesh.github.io/) based, software for continuum simulations. The package is designed to be easy to use, extendable and scalable. It currently includes two modules, [**COPSS-Hydrodynamics**](https://bitbucket.org/COPSS/copss-hydrodynamics-public.git) to solve hydrodynamic interactions in colloidal suspensions and [**COPSS-Polarization**](https://bitbucket.org/COPSS/copss-polarization-public) to solve electrostatic interactions between dielectric particles. The algorithms beneath **COPSS** have been published or under-review, but the code framework, user-interface, etc., are still rough. We are working on improving **COPSS** and appreciate your contributions.  


**COPSS-Hydrodynamics**
----------------------------

**COPSS-Hydrodynamics** solves the hydrodynamic interactions in colloidal suspensions by directly solve the Stokes flow.  It is based on an efficient $O(N)$ computational approach to model the dynamics of hydrodyna- -mically interacting Brownian or micron-sized particles in arbitrary geometries. A parallel finite element Stokes' solver is the center of the algorithm. 


**Installation**
-------------------------------------------
**COPSS** is written on [LIBMESH](http://libmesh.github.io/) framework. It also requires [PETSc](https://www.mcs.anl.gov/petsc/index.html) for parallel linear equation solvers and [SLEPc](http://slepc.upv.es/) for scalable Eigenvalue computations. Before installing **COPSS**, you need to install **LIBMESH** together with **PETSc** and **SLEPc**. To achieve the best parallel performance of COPSS, we suggest install it on a Linux cluster environment.
###0. System environment prep
Load or compile

 - [CMAKE](https://cmake.org/) (e.g., but not necessarily, version 3.6.2)
 - [GCC](https://gcc.gnu.org/) (e.g., but not necessarily, version 6.2)
 - [PYTHON](https://www.python.org/) (python 2)
 - [OPENMPI](https://www.open-mpi.org/) (e.g., but not necessarily, version 2.0.1)

###1. Install PETSc
 - Download PETSC's latest release ( version 3.7.4 or later ) from [PETSc download](https://www.mcs.anl.gov/petsc/download/index.html), or git clone PETSc repository:
	 
	 `mkdir $HOME/projects/`
	 
	 `cd $HOME/projects/`
	 
	 `git clone -b maint https://bitbucket.org/petsc/petsc petsc`
 - Configure PETSc: 
 
	 `cd $HOME/projects/petsc`
	 
	 `./configure --with-cc=mpicc --with-cxx=mpicxx --with-mpiexec=mpiexec --with-fc=mpif90 –download-fblaslapack --download-scalapack --download-mumps --download-superlu_dist --download-hypre --download-ml --download-parmetis --download-metis --download-triangle --download-chaco --with-debugging=0`
	 
	 **And then follow the instructions on screen to install and test the package.** 
 - Export environment variables:
	
	   `export PETSC_DIR=*/path/to/PETSC*`
	   `export PETSC_ARCH=*PETSC_ARCH_NAME*`
	
	  **Add the above codes to ~/.bashrc and `source ~/.bashrc` before next step. (`*/path/to/PETSC*` and `*PETSC_ARCH_NAME*` can be found on the screen after installation.)**

	If you meet any trouble, please refer to [PETSC installation](https://www.mcs.anl.gov/petsc/documentation/installation.html).
	
###2. Install SLEPc
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

###3. Install LIBMESH
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

	If you meet any trouble, please refer to [LIBMESH installation](https://libmesh.github.io/installation.html), or reach out to **LIBMESH** community for help.

###4. Install COPSS-hydrodynamics

 - Download the latest COPSS codes
 
 
	 `cd /path/to/where/you/want/to/install/copss`


	 `git clone https://bitbucket.org/COPSS/copss-hydrodynamics-public.git`



 __ ------------------On Master branch (only point particle systems are available)---------------------------------------__


 - Compile the codes (on Master branch)	 


	 `cd /path/to/copss/src/` 

	 `make`



 - Run an diffusion example of beads


 	  `cd /path/to/copss/examples/general_point_particle/beads/`

      `cp /path/to/copss/src/example-opt $PWD`	 

	  `cp /path/to/copss/tools/run.sh $PWD` 

      `bash run.sh` (You can define how many cores you want to run on in **run.sh**)
	  
	
	You need to set up your system in **point_particle_control.in** and **datafile (e.g., point_particle_data.in)**. More details can be found in our documentation.


 __ ------------------On Development branch (both PointParticle and RigidParticle are avaiable)---------------------------------------__

 - Compile the codes
 
	 `cd /path/to/copss/src/` 
	 
	 
	 `make package=POINTPARICLE (for point particle systems)  Or make package=RIGIDPARTICLE (for rigid particle systems)`	
	 
	 

 - Run an sedimentation example of rigid particles
 
 
 	  `cd /path/to/copss/examples/general_rigid_particle/sedimentation_benchmark/`
	  

       `modify run.sh dependending on your system`


      `bash run.sh`
	  

**Build documentation**
-------------------------------------------
After you have build **COPSS-Hydrodynamics** you can further build the documentation from doxygen/ directory. Make sure you have [Doxygen](http://www.stack.nl/~dimitri/doxygen/) ready:

    doxygen Doxyfile.bak

Then you can view the documentation in IE browser:

    google-chrome html.index.html

**What's next**
-------------------------------------------

 - Example systems for **finite size particles** using **Immersed Boundary Method (IBM)**.
 - Coupling COPSS with advanced sampling package [SSAGES](https://github.com/MICCoM/SSAGES-public) to enable advanced sampling for stokes systems ( **super exciting about this!!!** )
 - Coupling COPSS with dielectric polarization calculations ([COPSS-polarization](https://bitbucket.org/COPSS/copss-polarization-public) or [Image method](http://www.sciencedirect.com/science/article/pii/S0021979716301138) developed by COPSS team).
 - More user friendly interface. 
 - etc...

**Contribution**
-------------------------------------------
Hydrodynamic and electrostatic interactions are ubiquitous in nature and technology, including but not limited to colloidal, biological, granular systems, etc. **COPSS** is trying to solve these problems in a efficient and scalable manner. Contributions on the code structure, parallel-performance, applications are much appreciated. 

**Contact**
-------------------------------------------
 
 -  If you need help, have any questions or comments, join our mailing list [copss-users@googlegroups.com](Link URL)

     **GMail users**: just click "Join group" button

     **Everyone else**: send an email to [copss-users+subscribe@googlegroups.com](Link URL)

**Algorithm**
-------------------------------------------
 Our method paper has details on the implementation and tests of this code:  **Parallel $O(N)$ Stokes' solver towards scalable Brownian dynamics  in general geometries, submitted, 2017**


**Main contributors**
-------------------------------------------
 

 - [**Jiyuan Li**](https://scholar.google.com/citations?user=XE6JtJwAAAAJ&hl=en), Institute for Molecular Engineering, the University of Chicago. [LinkedIn](https://www.linkedin.com/in/jyliuchicago/) , [EMAIL](jyli@uchicago.edu)
 
 - [**Dr. Xikai Jiang**](https://www.researchgate.net/profile/Xikai_Jiang), Institute for Molecular Engineering, the University of Chicago. [EMAIL](xikai@uchicago.edu)
 
 - **Dr. Xujun Zhao**
 

**License**
-------------------------------------------
 
* The codes are open-source and distributed under the GNU GPL license, and may not be used for any commercial or for-profit purposes without our permission.
