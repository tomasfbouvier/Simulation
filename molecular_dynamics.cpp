#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include "molecular_dynamics.h"


double E=T*kB*N_ATOMS;

std::vector<Atom> generate_fcc_lattice(int N)
{
    std::vector<Atom> atoms;

    // Calculate number of atoms per side of lattice
    int n_per_side = static_cast<int>(std::cbrt(N));

    // Calculate lattice constant
    double lattice_constant = box_size / n_per_side;

    // Generate FCC lattice with N atoms
    for (int i = 0; i < n_per_side; i++)
    {
        for (int j = 0; j < n_per_side; j++)
        {
            for (int k = 0; k < n_per_side; k++)
            {
                Atom atom;
                atom.symbol = "Si";
                atom.x = i * lattice_constant;
                atom.y = j * lattice_constant;
                atom.z = k * lattice_constant;
                
                atoms.push_back(atom);
            }
        }
    }

    // Initialize velocities to set net force and momentum to zero
    std::random_device rd;
    std::mt19937 rng(rd());
    
    double sum_vx = 0, sum_vy = 0, sum_vz = 0;
    for (Atom& atom : atoms)
    {
        
        //std::cout<<"MyASDASDASD"<<std::sqrt(kB*T/mass)<<"\n";
        std::normal_distribution<double> dist_vx(0, std::sqrt(kB*T/mass));
        std::normal_distribution<double> dist_vy(0, std::sqrt(kB*T/mass));
        std::normal_distribution<double> dist_vz(0, std::sqrt(kB*T/mass));

        atom.vx = dist_vx(rng);
        atom.vy = dist_vy(rng);
        atom.vz = dist_vz(rng);

        sum_vx += atom.vx;
        sum_vy += atom.vy;
        sum_vz += atom.vz;
    }
    
    correct_energies(atoms);

    for (Atom& atom : atoms)
    {
        atom.vx -= sum_vx / N;
        atom.vy -= sum_vy / N;
        atom.vz -= sum_vz / N;
    }
    
    calculate_forces(atoms);
    
    for (Atom& atom : atoms)
    {
        atom.ax_prev = atom.ax;
        atom.ay_prev = atom.ay;
        atom.az_prev = atom.az;
    }
    
    return atoms;
}

void calculate_forces(std::vector<Atom>& atoms)
{
// Calculate forces acting on each atom
    for (Atom& atom : atoms)
    {
        atom.ax = 0;
        atom.ay = 0;
        atom.az = 0;
    }

    for (int i = 0; i < N_ATOMS -1; i++)
    {
        for (int j = i + 1; j < N_ATOMS; j++)
        {
            
            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            double dz = atoms[j].z - atoms[i].z;
            
            // Account for periodic boundary conditions
            
            if (dx >= box_size/2.) dx -= box_size;
            if (dx < -box_size/2.) dx += box_size;
            if (dy >= box_size/2.) dy -= box_size;
            if (dy < -box_size/2.) dy += box_size;
            if (dz >= box_size/2.) dz -= box_size;
            if (dz < -box_size/2.) dz += box_size;
            
            double r2 = dx*dx + dy*dy + dz*dz;
            double F =  24 * epsilon/(mass*sigma*sigma) * (2 * std::pow(sigma*sigma/r2, 7) - std::pow(sigma*sigma/r2, 4));
            
            atoms[i].ax -= F * dx;
            atoms[i].ay -= F * dy;
            atoms[i].az -= F * dz;
            
            atoms[j].ax += F * dx;
            atoms[j].ay += F * dy;
            atoms[j].az += F * dz;

        }
    }
}


void update_atoms_verlet(std::vector<Atom>& atoms)
{
    // Calculate forces acting on each atom
    // Update positions and velocities using velocity Verlet integration
    for (Atom& atom : atoms)
    {
        // Update positions
        atom.x += atom.vx * dt + 0.5 * atom.ax * dt * dt;
        atom.y += atom.vy * dt + 0.5 * atom.ay * dt * dt;
        atom.z += atom.vz * dt + 0.5 * atom.az * dt * dt;
    }

    // Keep inside the box
    apply_periodic_boundary_conditions(atoms);
    
    // Calculate new forces
    calculate_forces(atoms);
    
    for (Atom& atom : atoms)
    {
        // Update velocities
        atom.vx += 0.5 * (atom.ax + atom.ax_prev) * dt;
        atom.vy += 0.5 * (atom.ay + atom.ay_prev) * dt;
        atom.vz += 0.5 * (atom.az + atom.az_prev) * dt;
        
        atom.ax_prev = atom.ax;
        atom.ay_prev = atom.ay;
        atom.az_prev = atom.az;
    }
}

void update_atoms_euler(std::vector<Atom>& atoms)
{
    for (Atom& atom : atoms)
    {
        // Update velocity
        atom.vx += atom.ax * dt;
        atom.vy += atom.ay * dt;
        atom.vz += atom.az * dt;

        // Update position
        atom.x += atom.vx * dt;
        atom.y += atom.vy * dt;
        atom.z += atom.vz * dt;

        // Keep inside the box
        apply_periodic_boundary_conditions(atoms);
        
    }

}

double calculate_potential_energy(std::vector<Atom>& atoms)
{
    
    // Add correction term to potential energy
    double energy =  - N_ATOMS*(N_ATOMS-1) * 4. * M_PI * epsilon * sigma * sigma * sigma / (3 * r_cut * r_cut * r_cut);
    for (int i = 0; i < N_ATOMS -1; i++)
    {
        for (int j = i + 1; j < N_ATOMS; j++)
        {
            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            double dz = atoms[j].z - atoms[i].z;
            
            // Account for periodic boundary conditions
            
            if (dx >= box_size/2.) dx -= box_size;
            if (dx < -box_size/2.) dx += box_size;
            if (dy >= box_size/2.) dy -= box_size;
            if (dy < -box_size/2.) dy += box_size;
            if (dz >= box_size/2.) dz -= box_size;
            if (dz < -box_size/2.) dz += box_size;
            
            double r2 = dx*dx + dy*dy + dz*dz;
            
            if (r2 < r_cut*r_cut)
            {
                
                double r_6 = pow(sigma*sigma / r2, 3);
                double r_12 = r_6 * r_6;
                energy += 4 * epsilon * (r_12 - r_6);
            }
        }
    }
    return energy;
}


void correct_energies(std::vector<Atom>& atoms)
{
    double scale_factor=1.;
    kinetic_energy=0.;
    potential_energy=0.;
    for (int i = 0; i < N_ATOMS; i++)
    {
        kinetic_energy += 0.5 * mass * (atoms[i].vx * atoms[i].vx + atoms[i].vy * atoms[i].vy + atoms[i].vz * atoms[i].vz);
    }
    
    potential_energy = calculate_potential_energy(atoms);

    
    scale_factor = pow((E - potential_energy)/kinetic_energy , 0.5) ;
    
    for (int i = 0; i < N_ATOMS; i++)
    {
        atoms[i].vx *= scale_factor;
        atoms[i].vy *= scale_factor;
        atoms[i].vz *= scale_factor;
    }
    
    kinetic_energy*= scale_factor*scale_factor;
}


void apply_periodic_boundary_conditions(std::vector<Atom>& atoms)
{
    
    for (int i = 0; i < N_ATOMS; i++)
    {
        if (atoms[i].x < 0.) atoms[i].x += box_size;
        if (atoms[i].x > box_size) atoms[i].x -= box_size;
        if (atoms[i].y < 0.) atoms[i].y += box_size;
        if (atoms[i].y > box_size) atoms[i].y -= box_size;
        if (atoms[i].z < 0.) atoms[i].z += box_size;
        if (atoms[i].z > box_size) atoms[i].z -= box_size;
    }
}

void check_for_nan(const std::vector<Atom>& atoms)
{
    for (const Atom& atom : atoms)
    {
        if (std::isnan(atom.x) || std::isnan(atom.y) || std::isnan(atom.z))
        {
            std::cerr << "Error: Found nan value in atom positions!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (std::isnan(atom.vx) || std::isnan(atom.vy) || std::isnan(atom.vz))
        {
            std::cerr << "Error: Found nan value in atom velocities!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (std::isnan(atom.ax) || std::isnan(atom.ay) || std::isnan(atom.az))
        {
            std::cerr << "Error: Found nan value in atom accelerations!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

    }

    if (std::isnan(kinetic_energy) || std::isnan(potential_energy))
    {
        std::cerr << "Error: Found nan value in energy data!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


int main()
{
    // Set up constants


    // Generate FCC lattice
    std::vector<Atom> atoms = generate_fcc_lattice(N_ATOMS);

    // Open XYZ file for writing
    std::ofstream xyz_file("lattice.xyz");
    std::ofstream energies_file("energies.txt");
    // Evolve system for N_STEPS time steps
    for (int step = 0; step < N_STEPS; step++)
    {

        // Evolve system

        check_for_nan(atoms);
        update_atoms_verlet(atoms);
        
        std::cout<<step<<"\n";

        if (step % 100 == 0)
        {
            correct_energies(atoms);
            std::cout<<step<<"\t"<<"kinetic_energy: "<<kinetic_energy<<"\t"<<"potential_energy: "<<potential_energy<<"\t"<<"total_energy diff: "<<E-(kinetic_energy+potential_energy)<<"\n";
            // Do something with the energies (e.g. print them, store them, etc.)
        }

        // Write atoms to XYZ file
        xyz_file << N_ATOMS << std::endl;
        xyz_file << "Silicon FCC Lattice" << std::endl;
        for (const Atom& atom : atoms)
        {
            xyz_file << atom.symbol << " ";
            xyz_file << atom.x << " ";
            xyz_file << atom.y << " ";
            xyz_file << atom.z << std::endl;
        }
        
        
        kinetic_energy=0.;
        potential_energy=0.;
        for (int i = 0; i < N_ATOMS; i++)
        {
            kinetic_energy += 0.5 * mass * (atoms[i].vx * atoms[i].vx + atoms[i].vy * atoms[i].vy + atoms[i].vz * atoms[i].vz);

            
        }
        potential_energy = calculate_potential_energy(atoms);

        
        energies_file << kinetic_energy << "\t" << potential_energy << "\t"  << E << "\n";
    }
    
    

    // Close XYZ file
    xyz_file.close();

    return 0;
}
