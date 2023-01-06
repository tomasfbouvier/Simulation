# Molecular Dynamics Simulation
This is a molecular dynamics simulation in which a collection of atoms is simulated in a 3D space. The simulation starts with a FCC (face-centered cubic) lattice of atoms with a given number of atoms. The simulation then calculates the forces acting on each atom and updates the positions and velocities of the atoms over time using the velocity Verlet algorithm. The simulation also accounts for periodic boundary conditions, which means that if an atom moves outside of the boundaries of the 3D space, it will reappear on the opposite side as if the space were wrapped around like a torus. At each timestep, the forces on each atom are recalculated and the positions and velocities of the atoms are updated. The simulation continues for a given number of timesteps.

This code was written with the assistance of the artificial intelligence language model chatGPT.

## Lennard-Jones Potential

The Lennard-Jones potential is a model that describes the interactions between neutral atoms or molecules. It is given by the following formula:

$$ V_{LJ}(r) = 4 \epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right] $$

where $\epsilon$ is the depth of the potential well, $\sigma$ is the distance at which the potential energy is zero,
and $r$ is the distance between the two atoms. The first term represents the attractive interaction between the atoms,
while the second term represents the repulsive interaction. 
This potential is often used to model the interactions between atoms or molecules in materials science and statistical mechanics.

## Verlet algorithm 

The velocity Verlet algorithm is a popular method for integrating the equations of motion in molecular dynamics simulations. It allows the positions and velocities of the atoms to be updated in a stable and accurate way.

The velocity Verlet algorithm works by first updating the position of each atom using its current velocity and acceleration:

$$ \mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t) \Delta t + \frac{1}{2} \mathbf{a}(t) \Delta t^2 $$

where $\mathbf{r}(t)$ is the position of the atom at time $t$, $\mathbf{v}(t)$ is the velocity of the atom at time $t$, and $\mathbf{a}(t)$ is the acceleration of the atom at time $t$. $\Delta t$ is the time step of the simulation.

Next, the acceleration of each atom is updated using the forces acting on the atom:

$$ \mathbf{a}(t + \Delta t) = \frac{\mathbf{F}(t + \Delta t)}{m} $$

where $\mathbf{F}(t + \Delta t)$ is the force acting on the atom at time $t + \Delta t$, and $m$ is the mass of the atom.

Finally, the velocity of each atom is updated using the updated acceleration:

$$ \mathbf{v}(t + \Delta t) = \mathbf{v}(t) + \frac{1}{2} \left[ \mathbf{a}(t) + \mathbf{a}(t + \Delta t) \right] \Delta t $$

The above steps are then repeated for each time step of the simulation. This algorithm is relatively simple and efficient, making it a popular choice for molecular dynamics simulations.

## Example Movie of Molecular Dynamics in Action

![trial_movie](https://user-images.githubusercontent.com/57238320/211051267-f0db311e-78e2-427c-9617-73754a719b52.gif)



