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
    
    double kinv;
    double phivkinv;
    double phiv;
    double sqrdphivkinv;
    double phivv;
} box;


void verlet(box*, double, int);
double periodic(double, double);
double keep_inside_box(double , double );
void compute_acc(box**, int);
void correction(box*);


int main(int argc , char **argv){
	if(argc != 2){
		cout << " ---- Wrong number of arguments ---- " << endl;
		return -1;
	}
	const float n = atof(argv[1]);
    ifstream file;
	
	file.open("coordinates.txt", ios::in);
	
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



    ofstream energy_file;
	energy_file.open("results_energies.txt", ios::out);
	if(energy_file.fail()){
		cout<<"$$ | ENERGIES IS BROKEN | $$"<<endl;
		exit(1);
	}
	
 	ofstream meanval_file;
	meanval_file.open("mean_values.txt", ios::app);
	if(meanval_file.fail()){
		cout<<"$$ | MEANVAL IS BROKEN | $$"<<endl;
		exit(1);
	}	
	
	ofstream coordinates_file;
	coordinates_file.open("coordinates.txt", ios::out);
	if(coordinates_file.fail()){
		cout<<"Coordinates file is broken"<<endl;
	}
    
    
    
	double sumK=0;
    double sumU=0;
	double sumE=0;
	int a=0;
	        
    for (int t = 0 ; t < n ; t++){
        system.K = 0.0; 
        verlet(&system , 0.0001, t);
    /*    if(t==5000){
        	correction(&system);	
        }*/


        
///////////////////////record results after each 10 steps///////////////////////////        
        
        if(t%100==0){
         	sumK += system.K;
	   		sumU += system.U;
			energy_file << system.K << " " << system.U << " " << system.K + system.U << "\n";			
        }
        
    	if(t%1000==0){
    		cout<<"t= "<<t<<endl;
    		
    	}    

        	
        
    /*No me gusta nada este bucle pero bueno */
            

	   	a=t;
	    	
    }

	cout<<"tmax= "<<a<<endl;    
    
    meanval_file<<100*sumK/n<<" "<< 100*sumU/n<< " "<< 100*(sumK + sumU)/n<< " "<< 100*system.kinv/n<< " " <<100*system.phivkinv/n<< " " << 100*system.phiv/n<< " "<<100*system.sqrdphivkinv/n<< " "<<100*system.phivv/n<< "\n" ;
     
/* Record coordinates at the end to prepare next simulation */    

    for (int i = 0 ; i < N ; i++){
        coordinates_file<<system.world[i].r[0]<<"\t";
        coordinates_file<<system.world[i].r[1]<<"\t";
        coordinates_file<<system.world[i].r[2]<<"\t";
        coordinates_file<<system.world[i].p[0]<<"\t";
        coordinates_file<<system.world[i].p[1]<<"\t";
        coordinates_file<<system.world[i].p[2]<<"\t";
        coordinates_file<<system.world[i].a[0]<<"\t";
        coordinates_file<<system.world[i].a[1]<<"\t";
        coordinates_file<<system.world[i].a[2]<<"\n";
    }



/////////////Closes files/////////////
    
    energy_file.close();
    
    meanval_file.close();
    
    coordinates_file.close();

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

void verlet(box* BOX , double dt, int t){
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
    BOX->K=0; 
    for (int i = 0; i < N;i++){
        for (int k = 0; k<dim;k++){
            BOX->world[i].r[k] += BOX->world[i].p[k]*dt + 0.5*BOX->world[i].a[k]*dt*dt;
            BOX->world[i].r[k] = keep_inside_box(BOX->world[i].r[k] , BOX->length);
            BOX->world[i].p[k] += 0.5*BOX->world[i].a[k]*dt;
        }    
    }
    compute_acc(& BOX, t);
    for (int i = 0; i < N;i++){
        for (int k = 0; k<dim;k++){
            BOX->world[i].p[k] += 0.5*BOX->world[i].a[k]*dt;
        }
    }

    if(t%100==0){    
	    for (int i = 0 ; i < N ;i++){
	        for (int k = 0; k < dim ;k++){
				BOX->K += 0.5*(BOX->world[i].p[k])*(BOX->world[i].p[k]); 
	
		            
		    }
		}
  	}  
	///////////////////RECENTLY ADDED////////////////////////////    
    if(t%100==0){
	    BOX->kinv += 1/BOX->K;
	    BOX->phivkinv += BOX->phiv/BOX->K;
		BOX->sqrdphivkinv += BOX->phiv*BOX->phiv/BOX->K;
    }  

    
    
}

void compute_acc(box** BOX, int t){
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
    double U = - N*4*PI/(3.0*pow((*BOX)->length*0.5,3));
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

            for (int k = 0; k < dim; k++) {
                (*BOX)->world[i].a[k] += rij[k] * f;
                (*BOX)->world[j].a[k] -= rij[k] * f;
            }
            
            if(t%100==0){
               
	            (*BOX)->phiv += 1/(3*N/0.5)*f*pow(rsqd, 0.5);
	            
	            (*BOX)->phivv += 1/(3*N/0.5)*(1/(3*N/0.5)*(24 * (26 * pow(rsqd, -7) - 7*pow(rsqd, -4))*rsqd)-2*(*BOX)->phiv);
	            U += 4*(pow(rsqd, -6) - pow (rsqd, -3));
	       
            }
        }
        if(t%100==0){
        	(*BOX)->U = U;
        }

	}
}



