#include <stdio.h>      // Input and output
#include <stdlib.h>     // Pointers
#include <string.h>
#include <math.h>       // Math functions
#include <lapack.h>     // Linear algebra
#include "modelling_utils.h"
#define  PI 3.14159265358979323846264338328
#define dtrec 0.004
#define sec2ms 1.0e6
#define cfl 0.25
#define nderiv2 8


int main(){
       
// Declaracao de Variaveis
    
	// Parametros do modelo original 
    int nx,nz;    
	float dx,dz,x0,z0;
	float* vel;

	// Parametros do modelo extendido
	int nxx,nzz;
	struct b border;
	float* velextend;
	
	// Parametros de aquisição
	int it0;
	float xfmax,xfmin,dxf,zf,xg_offset_min,xg_offset_max,dxg,zg ;
	float freq,ttotal;
	float* ricker;

	// Parametros de estabilidade
	int ns,nt,ndt;
	float maxv,minv;
	float dxmax,dt;
	float deriv2[nderiv2]={-3.02359410, 1.75000000,-0.291666667, 0.0648148148,-0.0132575758,0.0021212122,-0.000226625227,0.0000118928690};
	float deriv2_sum=0;


	// Ler o modelo de velocidade
    FILE* velfile;
    FILE* velfile_hdr;
    FILE* velwrite;
    FILE* param;
	char* param_dir;
	char* velbin;
    char* velhdr;
    char* velextend_path;
    // Alocação de matrizes e vetores
    
    //Leitura do modelo de velocidade: hdr e bin
    velhdr = "/export/home/joan/ProjetoFWI/models/Vel.txt";
    velbin = "/export/home/joan/ProjetoFWI/models/Vel.bin";
    param_dir="/export/home/joan/ProjetoFWI/param.txt" ;
	velfile_hdr=fopen(velhdr,"r");    
    velfile=fopen(velbin,"rb");
	param=fopen(param_dir,"r");
	fscanf(param,"%f %f %f %f %f %f %f %f %f %f",&xfmax,&xfmin,&dxf,&zf,&xg_offset_min,&xg_offset_max,&dxg,&zg,&freq,&ttotal);
    printf("Print dos parametros:  %f \n %f \n %f \n %f \n %f \n %f \n %f \n %f \n  %f \n %f \n",
			xfmax,xfmin,dxf,zf,xg_offset_min,xg_offset_max,dxg,zg,freq,ttotal);
	fscanf(velfile_hdr,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx);
    
	vel=malloc(sizeof(float)*nx*nz);
	fread(vel,sizeof(float),nx*nz,velfile);
    fclose(velfile);
    fclose(velfile_hdr);
    //printf("nz= %i nx= %i z0=%f x0=%f dz=%f dx=%f \n",nz,nx,z0,x0,dz,dx);
    //printf("v0=%f v1=%f ", vel[1], vel[500]);   
  
    // Extending velocity model *******************************************
    border.nb=50;
	nxx=2*border.nb + nx;
    nzz=2*border.nb + nz;
    border.ixb=border.nb -1;
	border.izb=border.ixb;
	border.ize=nzz - border.nb;
   	
	velextend=malloc(sizeof(float)*nxx*nzz);
   	velextension(velextend,vel,nx,nz,border,nxx,nzz);
	
	velextend_path="./velextend.bin";
    velwrite=fopen(velextend_path,"wb");
    fwrite(velextend,4,nxx*nzz,velwrite);    
   


    //  Maximum and minimum of velocity model ******************************
	maxv=vel[0];
	minv=vel[0];
	for(int ii=0;ii<nx*nz;ii++){
		if(maxv < vel[ii]){
			maxv=vel[ii];
		}
		if(minv > vel[ii]){
			minv=vel[ii];
		}
	}

	// Setting parameters dt,dx for stability of FD

	ns= (int) (ttotal/dtrec) ;
	ns+=1;
	dxmax=minv/(6.0*freq);
	for(int ii=0;ii<nderiv2;ii++){
		deriv2_sum+=abs(deriv2[ii]);
	}
	dt=cfl*sqrt(2.0/deriv2_sum)*dxmax/maxv;
	ndt=(int)(dtrec/dt);
	ndt+=1;

	// Chosing dt for modelling to be a a divisor of dtrec
	if(ndt > 1){
		dt=(dtrec/ndt);
	}else{
		dt=dtrec;
	}


	nt=(int)(ttotal/dt);
	nt+=1;
	ricker=malloc(sizeof(float)*nt);
	
	it0=0;
	for(int it=it0;it<nt;it++){
		ricker[it]=fricker((it-it0)*dt,freq);
	}

	if(dx > dxmax){
		printf("Frequency is too high");
		return 0;
	}
	

	// Scaling  extended velocity for computational resource economy

	for(int ii=0;ii<nxx;ii++){
		for(int jj=0;jj<nzz;jj++){
			velextend[jj + ii*nzz]= pow(velextend[jj + ii*nzz],2.0)*pow(dt/dx,2.0);
		}
	}


	printf("Dt for stability: %.8f\n",dt);
	printf("Dx for stability:%f\n",dxmax);
	printf("ndt is:%i \n",ndt);
	printf("Maximum value of velocity: %f \n",maxv);
	printf("Minimum value of velocity: %f \n",minv);
	return 0;
	


}

