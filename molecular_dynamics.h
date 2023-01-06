//
//  molecular_dynamics.h
//  
//
//  Created by Tomás Fernández Bouvier on 2/1/23.
//

#ifndef generate_lattice_h
#define generate_lattice_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

const int N_ATOMS = 8000; // Number of atoms
const double mass = 39.948*1.66054e-27; // Mass of silicon atoms (Kg)
const double kB = 1.38064852e-23; // Boltzmann constant (J/K)
const double T = 1800; // Temperature (K)
const double dt = 1e-4; // Time step (ps)
const double box_size = 80.0 ; // Simulation box size (Å)
const int N_STEPS = 100; // Number of time steps to evolve system
const double r_cut = 100.; // (Å)

const double epsilon = 120 * kB;// 2.106e-20; // J
const double sigma = 3.40; // (Å)

double kinetic_energy;
double potential_energy;

struct Atom
{
std::string symbol;
double x, y, z;
double vx, vy, vz;
double ax, ay, az;
double ax_prev, ay_prev, az_prev;
};

std::vector<Atom> generate_fcc_lattice(int N);
void calculate_forces(std::vector<Atom>& atoms);
void correct_energies(std::vector<Atom>& atoms);
void update_atoms_verlet(std::vector<Atom>& atoms);
void update_atoms_euler(std::vector<Atom>& atoms);
double calculate_potential_energy(std::vector<Atom>& atoms);
void apply_periodic_boundary_conditions(std::vector<Atom>& atoms);
void check_for_nan(const std::vector<Atom>& atoms);

#endif /* generate_lattice_h */
