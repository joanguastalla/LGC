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

// Declaration of variables
    

// wavefield 
	float* p1,*p2;
// Parameters on original velocity model 
    int nx,nz;    
	float dx,dz,x0,z0;
	float* vel;

// Parameters on extended velocity model
	int nxx,nzz;
	struct b border;
	float* velextend;
	
// Aquisition parameters
	int it0;
	int isrc,isx,isz,idxs,nshots;
	int* igeo;
	float xsmax,xsmin,dxs,zs,xg_offset_min,xg_offset_max,dxg,zg ;
	float freq,ttotal;
	float* ricker;
	
// Stability parameters
	int ns,nt,ndt;
	float maxv,minv;
	float dxmax,dt;
	float deriv2[nderiv2]={-3.02359410, 1.75000000,-0.291666667, 0.0648148148,-0.0132575758,0.0021212122,-0.000226625227,0.0000118928690};
	float deriv2_sum=0;


// Directories for reading or writing files  
    char* param_dir;
	char* velbin;
    char* velhdr;
	char* rickerdir;
    char* velextend_path;
	param_dir="/export/home/joan/ProjetoFWI/param.txt" ;
    velbin = "/export/home/joan/ProjetoFWI/models/Vel.bin";
	velhdr = "/export/home/joan/ProjetoFWI/models/Vel.txt";
	velextend_path="./velextend.bin";
	rickerdir="/export/home/joan/ProjetoFWI/models/ricker.bin";
// File pointers to open files
	FILE* velfile;
    FILE* velfile_hdr;
    FILE* velwrite;
    FILE* param;
	FILE* ricker_bin;	
    // Alocação de matrizes e vetores
    
//Leitura do modelo de velocidade: hdr e bin
    velfile_hdr=fopen(velhdr,"r");    
    velfile=fopen(velbin,"rb");
	param=fopen(param_dir,"r");
	fscanf(param,"%f %f %f %f %f %f %f %f %f %f",&xsmax,&xsmin,&dxs,&zs,&xg_offset_min,&xg_offset_max,&dxg,&zg,&freq,&ttotal);
    printf("Print dos parametros:  %f \n %f \n %f \n %f \n %f \n %f \n %f \n %f \n  %f \n %f \n",
			xsmax,xsmin,dxs,zs,xg_offset_min,xg_offset_max,dxg,zg,freq,ttotal);
	fscanf(velfile_hdr,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx);
    
	vel=malloc(sizeof(float)*nx*nz);
	fread(vel,sizeof(float),nx*nz,velfile);
    fclose(velfile);
    fclose(velfile_hdr);
  
// Extending velocity model *******************************************
    border.nb=50;
	nxx=2*border.nb + nx;
    nzz=2*border.nb + nz;
    border.ixb=border.nb -1;
	border.izb=border.ixb;
	border.ize=nzz - border.nb;
   	
	velextend=malloc(sizeof(float)*nxx*nzz);
   	velextension(velextend,vel,nx,nz,border,nxx,nzz);
	
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
	
	// Half duration of ricker number of indexes it0	
	it0 = floor(sqrt(6.0/PI)/(freq*dt)); 
	for(int it=0;it<nt;it++){
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
// Wavefield temporal evolution common-shot situation
	p1=calloc(nxx*nzz,sizeof(float));
	p2=calloc(nxx*nzz,sizeof(float));

    ricker_bin=fopen(rickerdir,"wb");
    fwrite(ricker,sizeof(float),nt,ricker_bin);
	printf("nt:%d",nt);

// Index of source in the model
	
	isx=(int) (xsmin/dx);   // starting at 0 sample
	isz=(int) (zs/dz);	   // starting at 0 sample
	isrc= isz + (border.izb+1) + nzz*(isx + border.ixb +1);
	idxs=(int) (dxs/dx);
	idxs=idxs*nzz;
	nshots= (int) ((xsmax - xsmin)/dxs);
	nshots+=1;
	for(int is=0;is<nshots;is++){
		isrc=isrc + idxs*is;
		for(int it=0;it<nt;it++){
				p2[isrc]+= pow(dx,2)*ricker[it]*velextend[isrc];
	// Laplacian and field update
		}	
	}
	
	return 0;
	

}

