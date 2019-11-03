#include<iostream>
#include<fstream>
#include<string.h>
#include <stdlib.h>
#include<math.h>
#define N 500
#define PI 3.14159265359
#define E -605
#define dim 3
using namespace std;

typedef struct Particle{
    double r[dim];
    double p[dim];
    double a[dim];
} particle;

typedef struct Box{
    particle world[N];
    double length;
    double K;
    double U;
} box;


void verlet(box*, double);
double periodic(double, double);
double keep_inside_box(double , double );
void compute_acc(box**);
void correction(box*);


int main(int argc , char **argv){
	if(argc != 2){
		cout << " ---- Wrong number of arguments ---- " << endl;
		return -1;
	}
	const int n = atoi(argv[1]);
    cout << "  . \n ..^____/ \n`-. ___ ) \n ||  || VIVA FRANCO | PLVS VLTRA" << endl;
	ifstream file;
	
	file.open("results_cor.txt", ios::in);
	
	if(file.fail()){
		cout<< "/$$ -- FILE NOT FOUND -- $$/ " <<endl;
		exit(1);
	}
    box system;
    system.length = 10.0;
    for (int i = 0 ; i < N ; i++){
        file>>system.world[i].r[0];
        file>>system.world[i].r[1];
        file>>system.world[i].r[2];
        file>>system.world[i].p[0];
        file>>system.world[i].p[1];
        file>>system.world[i].p[2];
        file>>system.world[i].a[0];
        file>>system.world[i].a[1];
        file>>system.world[i].a[2];
    }
    file.close();
    double energies[n][3];
    for (int t = 0 ; t < n ; t++){
        system.K = 0.0;
        verlet(&system , 0.0001);
        if(t==150000){
        	correction(&system);	
        }
        energies[t][0] = system.K;
        energies[t][1] = system.U;
        energies[t][2] = system.K + system.U;
        cout << "Iteration $$$ ::: --->> " << t << endl;
        cout << "E -> " << system.U + system.K<< endl;
        cout << "K -> " << system.K << endl;
        cout << "U -> " << system.U << endl;
    }
    
    ofstream energy_file;
	energy_file.open("results_energies.txt", ios::out);
	if(energy_file.fail()){
		cout<<"$$ | FILE IS BROKEN | $$"<<endl;
		exit(1);
	}

    for (int t = 0 ; t < n ; t++){
        energy_file<<energies[t][0] << " " <<energies[t][1] << " " <<energies[t][2] << "\n";
    }
    energy_file.close();
    return 0;
}

double periodic(double u , double L){
    /*
    ---------------------------------
    -> Description:
    -> Arguments:
    @u double: Distance between two particles
    @L double: Length of the box
    -> Output:
    - double: Maps the distance such that it preserves systems periodicity 
    ---------------------------------
     */
    while (u < -0.5*L){
        u += L;
    }
    while (u >= 0.5*L){
        u -= L;
    }
    return u;
}

double keep_inside_box(double u , double L){
    /*
    ---------------------------------
    -> Description:
    -> Arguments:
    @u double: Position of a particle
    @L double: Length of the box
    -> Output:
    - double: Maps the distance such that it prevents a particle from escaping from the Box 
    ---------------------------------
     */
    if (u < 0.0){
        u += L;
    }
    else if (u > L){
        u -= L;
    }
    return u;
}

void correction(box* BOX ){
    /*
    ---------------------------------
    -> Description: Corrects the system after t steps in order to prevent the bounce effect in the total energy
    -> Arguments:
    @BOX box*: 
    -> Output:
    - void:
    ---------------------------------
     */	
	
	double KQ= E- BOX->U;
	double d=KQ/BOX->K;

	for(int i=0; i<N; i++){
		for(int k=0; k<3; k++){
			BOX->world[i].p[k] *= pow(d, 0.5);
		}
	}
	return;
}

void verlet(box* BOX , double dt){
    /*
    ---------------------------------
    -> Description: Integrates newton's equations of motion 
    -> Arguments:
    @box* BOX: pointer to the struct containing all the system
    @dt double: timestep
    -> Output:
    - void
    ---------------------------------
     */
    for (int i = 0; i < N;i++){
        for (int k = 0; k<dim;k++){
            BOX->world[i].r[k] += BOX->world[i].p[k]*dt + 0.5*BOX->world[i].a[k]*dt*dt;
            BOX->world[i].r[k] = keep_inside_box(BOX->world[i].r[k] , BOX->length);
            BOX->world[i].p[k] += 0.5*BOX->world[i].a[k]*dt;
        }    
    }
    compute_acc(& BOX);
    for (int i = 0; i < N;i++){
        for (int k = 0; k<dim;k++){
            BOX->world[i].p[k] += 0.5*BOX->world[i].a[k]*dt;
        }
    }
    for (int i = 0 ; i < N ;i++){
        for (int k = 0; k < dim ;k++){
            BOX->K += 0.5*(BOX->world[i].p[k])*(BOX->world[i].p[k]);
        }
    }
}

void compute_acc(box** BOX){
    /*
    ---------------------------------
    -> Description: Calculates forces among particles within the system 
    -> Arguments:
    @BOX box**: pointer to the pointer containing the struct to enable members to be updated
    -> Output:
    - void
    ---------------------------------
     */
    double rij[dim];
    double rsqd,f;
    double U = - N*4*PI/(3.0*pow((*BOX)->length/2.0,3));
    for (int i = 0; i < N; i++) {  
        for (int k = 0; k < dim; k++) {
            (*BOX)->world[i].a[k] = 0;
        }
    }
    for (int i = 0 ;i<N-1; i++){
        for (int j = i+1 ; j<N;j++){
            rsqd=0;
            for (int k = 0;k < dim;k++){
                rij[k] = periodic((*BOX)->world[i].r[k] - (*BOX)->world[j].r[k] , (*BOX)->length);
                rsqd += rij[k]*rij[k]; 
                }
            f = 24 * (2 * pow(rsqd, -7) - pow(rsqd, -4));
            U += 4*(pow(rsqd, -6) - pow (rsqd, -3));
            for (int k = 0; k < dim; k++) {
                (*BOX)->world[i].a[k] += rij[k] * f;
                (*BOX)->world[j].a[k] -= rij[k] * f;
                }
            }
        }
    (*BOX)->U = U;
}

