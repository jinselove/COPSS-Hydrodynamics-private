**COPSS: Force field**
=====================
1. Particle - particle force types
-------------------------------------------
 Particle-particle force types defines the force types between particle and particles. 

>     particle_particle_force_types = 'pp_ev_gaussian, pp_ev_gaussian_polymerChain, ...'
>     pp_ev_gaussian = 'param1, param2, ...'
>     pp_ev_gaussian_polymerChain = 'param1, param2, ...'

Supposing that we have two particles $i$ and $j$, located at $R_i$ and $R_j$ and the forces on which are $f_i$ and $f_j$ respectively. 

$\vec{f}_{ij}$: force acting on particle $i$ by particle $j$.
$\vec{R}_{ij}$: vector pointing from $i$ to $j$ , i.e., $\vec{R}_{ij} = \vec{R}_j - \vec{R}_i$, which is automatically updated due to periodic boundary conditions.
$\vec{r}_{ij}$ : unit vector of $\vec{R}_{ij}$
$R_{size}$: length of $\vec{R}_{ij}$, i.e., $\vec{R}_{ij} = R_{size}*\vec{r}_{ij}$
$a$: bead radius. All lengths are non-dimensionalized by this length. 
$b_k$: Kuhn length
$N_{k,s}$: number of Kuhn length per spring
$q_0$: maximum spring length, $q_0 = N_{k,s} * b_k$
$L$: contour length of the DNA molecule, $L = N_s * q_0$
$S_s^2$: radius of gyration of an ideal chain consisting of $N_{k,s} $ Kuhn segments, $S_s^2 = N_{k,s}*b_k^2/6$ 


----------
**pp_ev_gaussian**

> pp_ev_gaussian = '$c_1$, $c_2$'

pp_ev_gaussian defines a gaussian potential between point particles (**beads only**), two nondimensional parameters need to be given for this force type, $c_1$ (energy) and $c_2$ (length).

 Then this gaussian force:
> $\vec{f}_{ij} = -c_1*c_2*e^{-c_2*R_{size}^2}*\vec{r}_{ij} $
> $\vec{f}_i  += \vec{f}_{ij}$


----------
**pp_ev_gaussian_polymerChain**

> pp_ev_gaussian_polymerChain = '$ev$'

pp_ev_gaussian_polymerChain defines a gaussian potential between beads of worm-like polymer chain **(polymer chain only)**, the only required parameter $ev$ is the nondimensional excluded volume of beads. The coefficient of this gaussian potential is set by default as:

> $c_1 = ev * a ^3 * N_{k,s}^2 * (\frac{3.}{4. * \pi*S_s^2})^{3/2}$
> $c_2 =  3. * \frac{a^2}{4. * S_s^2}$

 Then this gaussian force:
> $\vec{f}_{ij} = -c_1*c_2*e^{-c_2*R_{size}^2}*\vec{r}_{ij} $
> $\vec{f}_i  += \vec{f}_{ij}$
	


----------
**pp_ev_lj_cut**

> pp_ev_lj_cut = '$\epsilon$, $\sigma$, $r_{cut}$'

pp_ev_lj_cut defines a Lennard-Jones potential between two particle $i$ and $j$. Three non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle diameter or slighter bigger, e.g., 2.1), $r_{cut}$ (cutoff radius) are required for this force field. 

Then the lj force:

>if  $R_{size} <=  r_{cut}$:
> $\vec{f}_{ij} = -24 * \epsilon * (2*(\frac{\sigma}{r_{ij}})^{12} - (\frac{\sigma}{r_{ij}})^{6} )$
> $\vec{f}_i  += \vec{f}_{ij}$
> else:
> $\vec{f}_i  += \vec{0}$

**pp_ev_lj_repulsive**

> pp_ev_lj_repulsive = '$\epsilon$, $\sigma$'

pp_ev_lj_repulsive defines a repulsive Lennard-Jones potential between two particle $i$ and $j$. Two non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle diameter or slighter bigger, e.g., 2.1) are required for this force field. 

$r_{cut}$ is set to be the equilibrium length where lj force is zero: 
> $r_{cut} = 2^{\frac{1.}{6.}} * \sigma$

Then the repulsive lj force:
>if  $R_{size} <=  r_{cut}$:
> $\vec{f}_{ij} = -24 * \epsilon * (2*(\frac{\sigma}{r_{ij}})^{12} - (\frac{\sigma}{r_{ij}})^{6} )$
> $\vec{f}_i  += \vec{f}_{ij}$
> else:
> $\vec{f}_i  += \vec{0}$


----------


**pp_ev_harmonic_repulsive**

> pp_ev_harmonic_repulsive = '$k$, $r_0$'

pp_ev_harmonic_repulsive defined a repulsive harmonic potential between particle $i$ and $j$. Two non-dimensional parameters, $k$(energy) and $r_0$ (equilibrium length) are required for this force field. 
Then the repulsive harmonic force:

> if $R_{size} < r_0$ : 
> $\vec{f}_{ij} = k * (R_{size} - r_0) * \vec{r}_{ij}$
> $\vec{f}_i  += \vec{f}_{ij}$
> else : 
> $\vec{f}_i  += \vec{0}$


----------


**pp_wormLike_spring**

> pp_wormLike_spring

pp_wormLike_spring defines spring forces for worm-like bead spring chains (**polymer chain only**). All parameters are set by default in COPSS.

> $c_1 = \frac{a}{2* b_k} $
> $L_s = \frac{N_{k,s}*b_k}{a}$ 

Then the spring force:

> $\vec{f}_{ij} = c_1*((1-\frac{R_{size}}{L_s})^{-2}-1.+4*\frac{R_{size}}{Ls})*\vec{r}_{ij}$
> $= \frac{a}{2*b_k}((1-\frac{R_{size}}{N_{k,s}*b_k/a})^{-2}-1.+4*\frac{R_{size}}{N_{k,s}*b_k/a}) *\vec{r}_{ij}$
> 
> $\vec{f}_i  += \vec{f}_{ij}$


----------


**p_constant**

> p_constant = '$f_x$, $f_y$, $f_z$'

p_constant defines a constant force field on all of the beads. Three parameters (force on $x, y, z$), $f_x, f_y, f_z$ are needed for the force field.
Then the constant force:

> $\vec{f}_{constant} = (f_x, f_y, f_z)$
> $\vec{f}_i += \vec{f}_{constant}$

 
 

2. Particle - wall force types
-----------------------------------------------------
 Particle-wall force types defines the force types between particles and wall, which has to be neither periodic boundary and inlet/outlet. 

>     particle_wall_force_types = 'pw_ev_empirical_polymerChain, pw_ev_lj_cut, ...'
>     pw_ev_empirical_polymerChain = 'param1, param2, ...'
>     pw_ev_lj_cut = 'param1, param2, ...'

Wall type can only be either **slit** or **sphere** for now, and will be extended to more types in further development. Supposing that we have particle $i$, located at $R_i$ and the forces on which is $f_i$.

$\vec{f}_{iw}$: force acting on particle $i$ by wall.
$\vec{R}_{iw}$: vector pointing from $i$ to wall.
> **if wall_type = 'slit'** :  $ \vec{R}_{i,lo} = \vec{box_{min} - \vec{R}_i},  \vec{R}_{i,hi} = \vec{box_{max} - \vec{R}_i} $ And we need to compute particle-wall interaction for lower wall and upper wall separately. 
> **if wall_type = 'sphere'** : $\vec{R}_{iw} = \vec{r}_i * (R_{sphere} - |\vec{R}_i|)$, where $\vec{r}_i$ is the unit vector of $\vec{R}_i$, $|\vec{R}_i|$ is the distance of particle $i$ to origin.  

$\vec{r}_{iw}$: unit vector of $\vec{R}_{iw}$.
$R_{size}$: length of $\vec{R}_{iw}$, i.e., $\vec{R}_{iw} = \vec{r}_{iw} * R_{size}$


----------


**pw_ev_empirical_polymerChain**

> pw_ev_empirical_polymerChain

pw_ev_empirical_polymerChain defines an empirical bead_wall repulsive potential on polymer beads (**polymer chain only**). All parameters are set by default in COPSS:

> $c_1 = a/b_k$
> $c2 = c1/\sqrt{N_{k,s}} = \frac{a}{b_k*\sqrt{N_{k,s}}}$
> $d_0 = 0.5/c_2 = \frac{b_k*\sqrt{N_{k,s}}}{2*a}$ 
> $c_0 = 25 * c_1 = \frac{25*a}{b_k}$

Then the empirical force:

> **if $R_{size} < d_0$**: 
> $\vec{f}_{iw} = -c_0 *(1- \frac{R_{size}}{d_0})^2*\vec{r}_{iw}$
> $= -\frac{25*a}{b_k}(1-\frac{2*R_{size}*a}{b_k*\sqrt{N_{k,s}}})^2*\vec{r}_{iw}$
> $\vec{f}_i += \vec{f}_{iw}$
> **else **:
> $\vec{f}_i += 0$

The corresponding potential is:

> **if $R_{size} < d_0$**:
> $U_i^{wall} = \frac{A_{wall}}{3*b_k/a * d_0}(R_{size} - d_0)^3$, where $A_{wall} = 25/a$
> **else**:
> $U_i^{wall} = 0$


----------


**pw_ev_lj_cut**

> pw_ev_lj_cut = '$\epsilon$, $\sigma$, $r_{cut}$'

pw_ev_lj_cut defines a Lennard-Jones potential between particle $i$ and the wall. Three non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle radius or slighter bigger, e.g., 1.05), $r_{cut}$ (cutoff radius) are required for this force field. 

Then the lj force:

>if  $R_{size} <=  r_{cut}$:
> $\vec{f}_{iw} = -24 * \epsilon * (2*(\frac{\sigma}{r_{iw}})^{12} - (\frac{\sigma}{r_{iw}})^{6} )$
> $\vec{f}_i  += \vec{f}_{iw}$
> else:
> $\vec{f}_i  += \vec{0}$


----------


**pw_ev_lj_repulsive**
> pw_ev_lj_repulsive = '$\epsilon$, $\sigma$'

pw_ev_lj_repulsive defines a repulsive Lennard-Jones potential between particle $i$ and the wall. Two non-dimensional parameters, $\epsilon$ (energy), $\sigma$ (particle radius or slighter bigger, e.g., 1.05) are required for this force field. 

$r_{cut}$ is set to be the equilibrium length where lj force is zero: 
> $r_{cut} = 2^{\frac{1.}{6.}} * \sigma$

Then the repulsive lj force:
>if  $R_{size} <=  r_{cut}$:
> $\vec{f}_{iw} = -24 * \epsilon * (2*(\frac{\sigma}{r_{iw}})^{12} - (\frac{\sigma}{r_{iw}})^{6} )$
> $\vec{f}_i  += \vec{f}_{iw}$
> else:
> $\vec{f}_i  += \vec{0}$


----------


**pw_ev_harmonic_repulsive**

 > pw_ev_harmonic_repulsive = '$k$, $r_0$'

pw_ev_harmonic_repulsive defined a repulsive harmonic potential between particle $i$ and the wall. Two non-dimensional parameters, $k$(energy) and $r_0$ (equilibrium length, e.g., 1.1) are required for this force field. 
Then the repulsive harmonic force:

> if $R_{size} < r_0$ : 
> $\vec{f}_{iw} = k * (R_{size} - r_0) * \vec{r}_{iw}$
> $\vec{f}_i  += \vec{f}_{iw}$
> else : 
> $\vec{f}_i  += \vec{0}$






