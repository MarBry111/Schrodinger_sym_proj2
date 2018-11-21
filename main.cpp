#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

FILE * file1;

int main() 
{

	file1 = fopen("aTextFile.txt","w");
	fclose(file1);

	int N = 100;
	double n =1;

	double kappa =0;
	double omega = 0;
	double tau = 1;

	double xd = 0.001;
	double x_k[N] = {};
	double psi_R[N] = {};
	double psi_I[N] = {};
	
	for(int i =0; i<=N; i++){
		x_k[i]=1./N*i;
		psi_R[i]=sqrt(2)*sin(n*M_PI*x_k[i]);
		psi_I[i]=0;
	}

	double H_R[N] = {};
	double H_I[N] = {};

	H_R[0] = 0;
	H_R[N] = 0;
	H_I[0] = 0;
	H_I[N] = 0;
	for(int i =1; i<=N-1; i++){
		H_R[i] = -0.5*(psi_R[i-1]+psi_R[i+1]-2*psi_R[i])/xd/xd + kappa*(x_k[i]-0.5)*psi_R[i]*sin(omega*tau);
		H_I[i] = -0.5*(psi_I[i-1]+psi_I[i+1]-2*psi_I[i])/xd/xd + kappa*(x_k[i]-0.5)*psi_I[i]*sin(omega*tau);
	}
	
	double dtau = 0.0001;
	int steps = 4000;

	for(int j = 0; j<= steps; j++){
		for(int i = 0; i<=N; i++){
			psi_R[i]=psi_R[i]+H_I[i]*dtau/2;
		}
		for(int i =1; i<=N-1; i++){
		H_R[i] = -0.5*(psi_R[i-1]+psi_R[i+1]-2*psi_R[i])/xd/xd + kappa*(x_k[i]-0.5)*psi_R[i]*sin(omega*tau);
		}
		for(int i = 0; i<=N; i++){
			psi_I[i]=psi_R[i]-H_R[i]*dtau;
		}
		for(int i =1; i<=N-1; i++){
		H_I[i] = -0.5*(psi_I[i-1]+psi_I[i+1]-2*psi_I[i])/xd/xd + kappa*(x_k[i]-0.5)*psi_I[i]*sin(omega*tau);
		}
		for(int i = 0; i<=N; i++){
			psi_R[i]=psi_R[i]+H_I[i]*dtau/2;
		}
		if(j%10==0){
			file1 = fopen("aTextFile.txt","a");
			fprintf(file1, "%f \n", H_I[1]);
			fclose(file1);
		}
	}


    return 0;
}

