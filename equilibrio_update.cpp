#include<iostream>
#include<fstream>
#include<string.h>
#include <stdlib.h>
#include<math.h>
#define N 500
#define PI 3.14159265359
#define E -605
using namespace std; 

typedef struct {
	double x[N];
	double y[N];
	double z[N];
	
	double vx[N];
	double vy[N];
	double vz[N];
	
	double ax_old[N];
	double ay_old[N];
	double az_old[N];

	double ax_new[N];
	double ay_new[N];
	double az_new[N];
	
	double pot[N];
	double Ec[N];
	
} Particulas;



void verlet(Particulas *par, ofstream archivo);
void potencial(Particulas ** par , int i , int j);


int main(int argc , char **argv){
	if(argc != 1){
		cout << " ---- Wrong number of arguments ---- " << endl;
	}
	const int n = atoi(argv[1]);

	ifstream archivo;
	
	archivo.open("results.txt", ios::in);
	
	if(archivo.fail()){
		cout<< "no se pudo abrir el archivo " <<endl;
		exit(1);
	}
	
	
	double aux[N][9];
	Particulas Part;
	//Esto puedo mejorarlo quitando el bucle aux
	for(int i=0; i<N; i++){
		for(int j=0; j<9; j++){
			archivo>>aux[i][j];
		}
		Part.x[i]= aux[i][0];
		Part.y[i]=aux[i][1];		
		Part.z[i]=aux[i][2];
		Part.vx[i]=aux[i][3];
		Part.vy[i]=aux[i][4];
		Part.vz[i]=aux[i][5];
		Part.ax_old[i]=aux[i][6];								
		Part.ay_old[i]=aux[i][7];
		Part.az_old[i]=aux[i][8];
		Part.pot[i] = 0.0;
	}
	archivo.close();
	
	ofstream archivo_energias;
	archivo_energias.open("results_energias.txt", ios::out);
	if(archivo_energias.fail()){
		cout<<"no se pudo abrir el archivo"<<endl;
		exit(1);
	}
	
	for(int i=0; i<n; i++){
		double  rc=5,	U = - N*4*PI/(3.0*pow(rc,3)) , K=0.0;
		
		verlet(&Part, archivo_energias);	

		//////////////// Corrección despues del bote /////////////////////	
		

///////////////////////ESTOS BUCLES SON HORRIBLES PERO NO QUEDA OTRA QUE TRABAJAR CON ELLOS PORQUE LA FUNCION NO ME LOS QUIERE DEVOLVER//////////////////////

		if(i==500){
			for(int l=0; l<N-1; l++){
				U+= Part.pot[l];
				
			}
			for(int l=0; l<N; l++){
				K+= Part.Ec[l];	
			}			
			
	
			
			double Tq=E-U;
			
			double d=Tq/K;
			
			
			cout<<"d= "<<d<<endl;
			for(int j=0; j<N; j++){
				
			(Part.vx[j]= Part.vx[j]*pow(d,0.5));
			(Part.vy[j]= Part.vy[j]*pow(d,0.5));
			(Part.vz[j]= Part.vz[j]*pow(d,0.5));
			
			
		/////////////////// YA NO PETA PERO NO ENTIENDO PQ NO DA !!!!!!!!!!!!/////////////////////////
		
		////////////////// FIN DEl INTENTO DE CORRECCIÓN///////////////////////////////////////////							
			}
			
		}

	}
	
	archivo_energias.close()


	ofstream archivo_prima;
	archivo_prima.open("results_prime.txt", ios::out);
	if(archivo_prima.fail()){
		cout<<"no se pudo abrir el archivo"<<endl;
		exit(1);
	}


	
	for(int i=0; i<N; i++){
		archivo_prima<<Part.x[i]<<"\t";
		archivo_prima<<Part.y[i]<<"\t";		
		archivo_prima<<Part.z[i]<<"\t";
		archivo_prima<<Part.vx[i]<<"\t";
		archivo_prima<<Part.vy[i]<<"\t";
		archivo_prima<<Part.vz[i]<<"\t";
		archivo_prima<<Part.ax_old[i]<<"\t";								
		archivo_prima<<Part.ay_old[i]<<"\t";
		archivo_prima<<Part.az_old[i]<<"\t";
		archivo_prima<<Part.pot[i]<<"\n";
	}
	archivo_prima.close();

		
	

	return 0;
}	
	
void verlet(Particulas * par, ofstream archivo){
///////////////// This function makes an input system evolve in time dt and returns both kinetic and potencial energy/////////// 
	
	double dt=0.0001;
	double rc= 5;
	double U, cin;
	double L=10;


 	U = - N*4*PI/(3.0*pow(rc,3));
    cin = 0.0;
	

	

	
	for(int i=0; i<N; i++){
	
		par->x[i] += (par->vx[i])*dt +0.5*(par->ax_old[i])*dt*dt;
		par->vx[i] += 0.5*(par->ax_old[i]+par->ax_new[i])*dt;
		par->y[i] += (par->vy[i])*dt +0.5*(par->ay_old[i])*dt*dt;
		par->vy[i] += 0.5*(par->ay_old[i]+par->ay_new[i])*dt;
		par->z[i] += (par->vz[i])*dt +0.5*(par->az_old[i])*dt*dt;
		par->vz[i] += 0.5*(par->az_old[i]+par->az_new[i])*dt;			
		par->ax_old[i] = par->ax_new[i];
		par->ay_old[i] = par->ay_new[i];
		par->az_old[i] = par->az_new[i];
		par->Ec[i] =  (par->vx[i]*par->vx[i] + par->vy[i]*par->vy[i]+par->vz[i]*par->vz[i])/2.0;
		cin += par->Ec[i];
		
		
		///////////////////paquete de instrucciones por si sale de la caja///////////////////
		if(par->x[i]<0){
			par->x[i] += L;
		}
		else{
			if(par->x[i]>L){
				par->x[i] -= L;
			}
				
		}
			
		if(par->y[i]<0){
			par->y[i] += L;
		}	
		else{
			if(par->y[i]>L){
				par->y[i] -= L;
					
			}
				
		}
			
			
		if(par->z[i]<0){
			par->z[i] += L;
		
		}
		else{
			if(par->z[i]>L){
				par->z[i] -= L;
					
			}
				
		}
				
		
	}
	
////////////////////// Este bucle no estoy seguro de que estea bien//////////////////
	
	for(int i=0; i<N-1; i++){
		par->pot[i]=0.0;
		for(int j=i+1; j<N; j++){
			potencial(&par , i , j);		
		}
		U += par->pot[i];
	}	
	


//////////////////////// Aquí pasa algo muy chungo /////////////////////////////
	
	archivo<<cin<<"\n";
	archivo<<U<<"\n";

	cout<<"E= "<<cin<<endl;
	


	return ;
}



void potencial(Particulas** par , int i , int j){
	double rc=5, L=10, Fmod=0;
	
	double rx= ((*par)->x[j] - (*par)->x[i]);
	double ry= ((*par)->y[j] - (*par)->y[i]);
	double rz= ((*par)->z[j] - (*par)->z[i]);			
			
	double alpha;
	if(rx<-rc){
		alpha=-1;
	}	
	else{
		if(rx<rc){
			alpha=0;
		}
		else{
			alpha= 1;
		}
	}	
	rx-= alpha*L;

	if(ry<-rc){
		alpha=-1;
	}	
	else{
		if(ry<rc){
			alpha=0;
		}
		else{
			alpha= 1;
		}
 								
	}
				
	ry-= alpha*L;		
	if(rz<-rc){
		alpha=-1;
	}	
	else{
		if(rz<rc){
			alpha=0;
		}
		else{
			alpha= 1;
		}
		
	}
				
	rz -= alpha*L;	
	
	double r= pow(rx*rx+ry*ry+rz*rz, 0.5);
	
	double aux= 4*(pow(1/r, 12) - pow (1/r, 6));
		
	(*par)->pot[i] += aux;
	(*par)->pot[j] += aux;
	

	Fmod = 24*(2*pow(r,-12)-pow(r, -6))*pow(r, -2);
	
/////////////////////// Empiezo a sospechar que aquí hay algo que quizás no está bien	
	
	(*par)->ax_new[i]+= Fmod*rx/r;
	(*par)->ay_new[i]+= Fmod*ry/r;
	(*par)->az_new[i]+= Fmod*rz/r;
	(*par)->ax_new[j]+= -Fmod*rx/r;
	(*par)->ay_new[j]+= -Fmod*ry/r;
	(*par)->az_new[j]+= -Fmod*rz/r;
	


	return;

}

	












