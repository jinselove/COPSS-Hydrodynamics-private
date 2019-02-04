.. _tutorials:

Tutorials
==========



Force Field
--------------

1. Particle - particle force types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 Particle-particle force types defines the force types between particle and particles:
 ::
 
    particle_particle_force_types = 'pp_ev_gaussian, pp_ev_gaussian_polymerChain, ...'
    pp_ev_gaussian = 'param1, param2, ...'
    pp_ev_gaussian_polymerChain = 'param1, param2, ...'


Supposing that we have two particles $i$ and $j$, located at $R_i$ and $R_j$ and the forces on which are $f_i$ and $f_j$ respectively.

$\vec{f}\ *{ij}$: force acting on particle $i$ by particle $j$.

$\vec{R}*\ {ij}$: vector pointing from $i$ to $j$ , i.e., $\vec{R}_{ij} = \vec{R}_j - \vec{R}\ *i$, which is automatically updated due to periodic boundary conditions.

$\vec{r}*\ {ij}$ : unit vector of $\vec{R}\ *{ij}$
$R*\ {size}$: length of $\vec{R}\ *{ij}$, i.e., $\vec{R}*\ {ij} = R\ *{size}*\vec{r}*\ {ij}$
$a$: bead radius. All lengths are non-dimensionalized by this length.
$b\ *k$: Kuhn length
$N*\ {k,s}$: number of Kuhn length per spring
$q_0$: maximum spring length, $q\ *0 = N*\ {k,s} * b_k$
$L$: contour length of the DNA molecule, $L = N_s * q_0$
$S\ *s^2$: radius of gyration of an ideal chain consisting of $N*\ {k,s} $ Kuhn segments, $S\ *s^2 = N*\ {k,s}*b_k^2/6$



* **pp_ev_gaussian**:
::

   pp_ev_gaussian = '$c_1$, $c_2$'


pp_ev_gaussian defines a gaussian potential between point particles (\ **beads only**\ ), two nondimensional parameters need to be given for this force type, $c_1$ (energy) and $c_2$ (length).
Then this gaussian force:

   $\vec{f}_{ij} = -c_1\ *c_2*\ e^{-c\ *2*R*\ {size}^2}*\vec{r}_{ij} $
   
   $\vec{f}\ *i  += \vec{f}*\ {ij}$


* **pp_ev_gaussian_polymerChain**:

::

   pp_ev_gaussian_polymerChain = '$ev$'


pp_ev_gaussian_polymerChain defines a gaussian potential between beads of worm-like polymer chain **(polymer chain only)**\ , the only required parameter $ev$ is the nondimensional excluded volume of beads.
The coefficient of this gaussian potential is set by default as:

   $c\ *1 = ev * a ^3 * N*\ {k,s}^2 * (\frac{3.}{4. * \pi\ *S_s^2})^{3/2}$
   $c_2 =  3. * \frac{a^2}{4. * S_s^2}$

Then this gaussian force:

   $\vec{f}_{ij} = -c_1\ *c_2*\ e^{-c\ *2*R*\ {size}^2}*\vec{r}_{ij} $
   $\vec{f}\ *i  += \vec{f}*\ {ij}$



* **pp_ev_lj_cut**:
::

   pp_ev_lj_cut = '$\epsilon$, $\sigma$, $r*\ {cut}$'


pp_ev_lj_cut defines a Lennard-Jones potential between two particle $i$ and $j$. Three non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle diameter or slighter bigger, e.g., 2.1), $r*\ {cut}$ (cutoff radius) are required for this force field.

Then the lj force:

   if  $R\ *{size} <=  r*\ {cut}$:
   $\vec{f}\ *{ij} = -24 * \epsilon * (2*(\frac{\sigma}{r*\ {ij}})^{12} - (\frac{\sigma}{r_{ij}})^{6} ) * r_{ij} / r_{ij}^2$

   $\vec{f}\ *i  += \vec{f}*\ {ij}$

   else:
   $\vec{f}_i  += \vec{0}$


* **pp_ev_lj_repulsive**:
::

   pp_ev_lj_repulsive = '$\epsilon$, $\sigma$'


pp_ev_lj_repulsive defines a repulsive Lennard-Jones potential between two particle $i$ and $j$. Two non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle diameter or slighter bigger, e.g., 2.1) are required for this force field.

$r_{cut}$ is set to be the equilibrium length where lj force is zero:

   $r_{cut} = 2^{\frac{1.}{6.}} * \sigma$


Then the repulsive lj force:

   if  $R\ *{size} <=  r*\ {cut}$:
   $\vec{f}\ *{ij} = -24 * \epsilon * (2*(\frac{\sigma}{r*\ {ij}})^{12} - (\frac{\sigma}{r_{ij}})^{6} ) * r_{ij} / r_{ij}^2$

   $\vec{f}\ *i  += \vec{f}*\ {ij}$

   else:
   $\vec{f}_i  += \vec{0}$


* **pp_ev_harmonic_repulsive**:
::

   pp_ev_harmonic_repulsive = '$k$, $r_0$'


pp_ev_harmonic_repulsive defined a repulsive harmonic potential between particle $i$ and $j$. Two non-dimensional parameters, $k$(energy) and $r_0$ (equilibrium length) are required for this force field.
Then the repulsive harmonic force:

   if $R_{size} < r\ *0$ :
   $\vec{f}*\ {ij} = k * (R_{size} - r_0) * \vec{r}_{ij}$
   $\vec{f}\ *i  += \vec{f}*\ {ij}$
   else :
   $\vec{f}_i  += \vec{0}$



* **pp_wormLike_spring**:
::

   pp_wormLike_spring


pp_wormLike_spring defines spring forces for worm-like bead spring chains (\ **polymer chain only**\ ). All parameters are set by default in COPSS.

   $c_1 = \frac{a}{2\ * b_k} $
   $L\ *s = \frac{N*\ {k,s}*\ b_k}{a}$


Then the spring force:

   $\vec{f}_{ij} = c\ *1*((1-\frac{R*\ {size}}{L\ *s})^{-2}-1.+4*\frac{R*\ {size}}{Ls})\ *\vec{r}_{ij}$
   $= \frac{a}{2*\ b\ *k}((1-\frac{R*\ {size}}{N_{k,s}\ *b_k/a})^{-2}-1.+4*\ \frac{R\ *{size}}{N*\ {k,s}\ *b_k/a}) *\ \vec{r}_{ij}$

   $\vec{f}\ *i  += \vec{f}*\ {ij}$


* **p_constant**:
::

   p_constant = '$f_x$, $f_y$, $f_z$'


p_constant defines a constant force field on all of the beads. Three parameters (force on $x, y, z$), $f_x, f_y, f_z$ are needed for the force field.
Then the constant force:

   $\vec{f}_{constant} = (f_x, f_y, f_z)$
   $\vec{f}\ *i += \vec{f}*\ {constant}$


2. Particle - wall force types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 Particle-wall force types defines the force types between particles and wall,
 which has to be neither periodic boundary and inlet/outlet:
 ::
 
      particle_wall_force_types = 'pw_ev_empirical_polymerChain, pw_ev_lj_cut, ...'
      pw_ev_empirical_polymerChain = 'param1, param2, ...'
      pw_ev_lj_cut = 'param1, param2, ...'


Wall type can only be either **slit** or **sphere** for now, and will be extended to more types in further development. Supposing that we have particle $i$, located at $R_i$ and the forces on which is $f_i$.

    $\vec{f}\ *{iw}$: force acting on particle $i$ by wall.
    $\vec{R}*\ {iw}$: vector pointing from $i$ to wall.


**if wall_type = 'slit'** :  $ \vec{R}\ *{i,lo} = \vec{box*\ {min} - \vec{R}\ *i},  \vec{R}*\ {i,hi} = \vec{box_{max} - \vec{R}_i} $ And we need to compute particle-wall interaction for lower wall and upper wall separately.
**if wall_type = 'sphere'** : $\vec{R}_{iw} = \vec{r}\ *i * (R*\ {sphere} - |\vec{R}_i|)$, where $\vec{r}_i$ is the unit vector of $\vec{R}_i$, $|\vec{R}_i|$ is the distance of particle $i$ to origin.


    $\vec{r}\ *{iw}$: unit vector of $\vec{R}*\ {iw}$.
    $R\ *{size}$: length of $\vec{R}*\ {iw}$, i.e., $\vec{R}\ *{iw} = \vec{r}*\ {iw} * R_{size}$


* **pw_ev_empirical_polymerChain**:
::

   pw_ev_empirical_polymerChain

pw_ev_empirical_polymerChain defines an empirical bead_wall repulsive potential on polymer beads (\ **polymer chain only**\ ). All parameters are set by default in COPSS:

   $c_1 = a/b\ *k$
   $c2 = c1/\sqrt{N*\ {k,s}} = \frac{a}{b\ *k*\sqrt{N*\ {k,s}}}$
   $d_0 = 0.5/c_2 = \frac{b\ *k*\sqrt{N*\ {k,s}}}{2\ *a}$
   $c_0 = 25 * c_1 = \frac{25*a}{b_k}$

Then the empirical force:

   **if $R_{size} < d_0$**\ :
   $\vec{f}_{iw} = -c\ *0 *(1- \frac{R*\ {size}}{d\ *0})^2*\vec{r}*\ {iw}$
   $= -\frac{25\ *a}{b_k}(1-\frac{2*\ R_{size}\ *a}{b_k*\ \sqrt{N\ *{k,s}}})^2*\vec{r}*\ {iw}$
   $\vec{f}\ *i += \vec{f}*\ {iw}$
   **else **\ :
   $\vec{f}_i += 0$


The corresponding potential is:

   **if $R_{size} < d_0$**\ :
   $U\ *i^{wall} = \frac{A*\ {wall}}{3\ *b_k/a * d\ *0}(R*\ {size} - d\ *0)^3$, where $A*\ {wall} = 25/a$
   **else**\ :
   $U_i^{wall} = 0$



* **pw_ev_lj_cut**:
::

   pw_ev_lj_cut = '$\epsilon$, $\sigma$, $r*\ {cut}$'


pw_ev_lj_cut defines a Lennard-Jones potential between particle $i$ and the wall. Three non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle radius or slighter bigger, e.g., 1.05), $r*\ {cut}$ (cutoff radius) are required for this force field.
Then the lj force:

   if  $R\ *{size} <=  r*\ {cut}$:
   $\vec{f}\ *{iw} = -24 * \epsilon * (2*(\frac{\sigma}{r*\ {iw}})^{12} - (\frac{\sigma}{r_{iw}})^{6} ) * r_{iw} / r_{iw}^2$

   $\vec{f}\ *i  += \vec{f}*\ {iw}$

   else:
   $\vec{f}_i  += \vec{0}$


* **pw_ev_lj_repulsive**:
::

   pw_ev_lj_repulsive = '$\epsilon$, $\sigma$'

pw_ev_lj_repulsive defines a repulsive Lennard-Jones potential between particle $i$ and the wall. Two non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle radius or slighter bigger, e.g., 1.05) are required for this force field.

$r_{cut}$ is set to be the equilibrium length where lj force is zero:

   $r_{cut} = 2^{\frac{1.}{6.}} * \sigma$

Then the repulsive lj force:

   if  $R\ *{size} <=  r*\ {cut}$:
   $\vec{f}\ *{iw} = -24 * \epsilon * (2*(\frac{\sigma}{r*\ {iw}})^{12} - (\frac{\sigma}{r_{iw}})^{6} ) * r_{iw} / r_{iw}^2$

   $\vec{f}\ *i  += \vec{f}*\ {iw}$

   else:
   $\vec{f}_i  += \vec{0}$


* **pw_ev_harmonic_repulsive**:
::

   pw_ev_harmonic_repulsive = '$k$, $r_0$'


pw_ev_harmonic_repulsive defined a repulsive harmonic potential between particle $i$ and the wall. Two non-dimensional parameters, $k$(energy) and $r_0$ (equilibrium length, e.g., 1.1) are required for this force field.
Then the repulsive harmonic force:

   if $R_{size} < r\ *0$ :
   $\vec{f}*\ {iw} = k * (R_{size} - r_0) * \vec{r}_{iw}$
   $\vec{f}\ *i  += \vec{f}*\ {iw}$
   else :
   $\vec{f}_i  += \vec{0}$
   
Integration tests
-------------------
The purpose of integration tests is to make sure new developments do not mess up the system. So far,
we have prepared several integration test systems:
    
    1) **PointParticle_Polymer_BD_HI**: Single polymer chain diffusing in a slit channel with HI considered.
    2) **RigidParticle_Sphere_Sedimentation_HI**: Single spherical particle sedimentate in a slit channel with HI considered.
    3) **PointParticle_Stokes_GGEM**: A validation system for Stokes GGEM. Three point forces are placed in\
    a infinite domain. The analytical solution at the box boundary is used as the boundary conditions. Solve
    the system using Stokes GGEM solver and compare the numerical velocity with the analytical velocity within\
    the domain.
     
The benchmark systems are located at **$COPSS_DIR/tests/integration_tests/resources/**. For each of the system, benchmark output
for integration tests were generated by running run.sh and input files in each of the folder. And the results are stored in the **/output**,
for example, "$COPSS_DIR/tests/integration_tests/resources/PointParticle_Polymer_BD_HI/output/".

The workflow in integration tests is:

    1) Check all documents, include "required_inputs" and "validation_outputs" defined in
    'test.json' are located in the respective resource folder
    2) The program enters respective resource folder for each test; 
    3) Compile the respective package; 
    4) Run the simulation using control file and data file; 
    5) Compare the simulation outputs and benchmark outputs in 'output/': 
    
        * If "cmp_method" == 'cmpfile', the program will compare the files directly.
        * If "cmp_method" == "cmpdf", the program will call the 'compare_output()' function
        defined in 'compare_output.py' is the respective test folder. One can customize the
        function for customized comparison of a test. 
        * Other "cmp_method" is not supported yet.
    6) Delete all intermediate files generated during tests and Leave the program 


How to run integration test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. change directory to integration test folder:
::
    
    cd $COPSS_DIR/tests/integration_tests
        
2. run the test script using Python 2.7.12 or up with at least 4 cpus available:
::
    
     python test.py
        
How to add a new integration test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For example, if we were to add a new test system called **NEW_TEST** to the integration test.

    1. Create a new folder:
    ::
    
        mkdir $COPSS_DIR/tests/integration_tests/resources/NEW_TEST
        
    2. Create corresponding **input files**, **"run.sh"** and **"zclean.sh"** under the new folder; If
    'cmp_method' = 'cmpdf', one needs to create a compare_output.py file in which a compare_output() function
    needs to be defined to customized file comparison.
    
    
    3. Run the simulation use the new files created in step 2 and store necessary outputs under folder, **$COPSS_DIR/tests/integration_tests/resources/NEW_TEST/output**
    
    4. Modify **$COPSS_DIR/tests/integration_tests/test.json** to include the new test.
    
    5. Testing the new integration test.

    
