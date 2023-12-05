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
#define nderiv2 3


int main(){

// Declaration of variables
    

// wavefield 
	float* p1,*p2;
	float* pfield;
	int iframe,im,ntrec;

// Parameters on original velocity model 
    int nx,nz;    
	float dx,dz,x0,z0;
	float* vel;

// Parameters on extended velocity model
	int nxx,nzz;
	struct b border;
	float* velextend;
	int i0_model,imodel;
	float laplacian,dswap;
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
	//float deriv2[nderiv2]={-3.02359410, 1.75000000,-0.291666667, 0.0648148148,-0.0132575758,0.0021212122,-0.000226625227,0.0000118928690};
	float deriv2[nderiv2]={-5/2,4/3,-1/12};
	float deriv2_sum=0;


// Directories for reading or writing files  
    char* param_dir;
	char* velbin;
    char* velhdr;
	char* rickerdir;
    char* velextend_path;
	char* moviedir;
	param_dir="/export/home/joan/ProjetoFWI/param.txt" ;
    velbin = "/export/home/joan/ProjetoFWI/models/Vel.bin";
	velhdr = "/export/home/joan/ProjetoFWI/models/Vel.txt";
	velextend_path="./velextend.bin";
	rickerdir="/export/home/joan/ProjetoFWI/models/ricker.bin";
	moviedir="/scratch/joan/wavemovie.bin";
// File pointers to open files
	FILE* velfile;
    FILE* velfile_hdr;
    FILE* velwrite;
    FILE* param;
	FILE* ricker_bin;	
    FILE* wavemovie;
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
  	for (int i = 0;i < nx*nz; ++i) {
  		vel[i]=2000;
  	}


// Extending velocity model *******************************************
    border.nb=150;
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
	dt=cfl*sqrt(2.0/deriv2_sum)*dx/maxv;
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
    i0_model=nzz*(border.ixb+1) + (border.izb + 1); 
	fwrite(ricker,sizeof(float),nt,ricker_bin);
	printf("nt:%d",nt);

// Index of source in the model
	
	isx=(int) (xsmin/dx);   // starting at 0 sample
	isz=(int) (zs/dz);	   // starting at 0 sample
	isrc= i0_model + isz +  nzz*isx;
	idxs=(int) (dxs/dx);
	idxs=idxs*nzz;
	nshots=(int) ((xsmax - xsmin)/dxs);
	nshots+=1;
// Wave movie variables 
wavemovie=fopen(moviedir,"wb");
ntrec=(int) (nt/ndt);
pfield=calloc(nx*nz,sizeof(float));
float modulo;
// Wave propagation in borders
int ix,iz;
float gamma_x[nxx];
float gamma_z[nzz];
float gamma, mgamma, invpgamma, beta;
beta=PI*freq*dt;
for(iz = 0; iz < nzz; ++iz) {
	gamma_z[iz]=0.0;
}
for (int ix = 0; ix < nxx; ++ix) {
	gamma_x[ix]=0.0;
}

for (ix = 0; ix < border.nb; ++ix) {
	gamma=beta*pow((((float) (ix+1))/((float) border.nb)),2.0);
	gamma_x[border.ixb - ix]=gamma;
	gamma_x[(border.ixb + nx + 1) + ix]=gamma;
}	

for (iz = 0; iz < border.nb; ++iz) {
	gamma=beta*pow((((float) (iz+1))/((float) border.nb)),2.0);
	gamma_z[border.izb - iz]=gamma;
	gamma_z[(border.izb + nz + 1) + iz]=gamma;
}	





// Wave evolution without borders
for(int is=0;is<nshots;is++){
		for(int iswap=0;iswap<nxx*nzz;iswap++){
					p1[iswap]=0.0;
					p2[iswap]=0.0;
		}
		isrc=isrc + idxs*is;
		iframe=0;
		for(int it=0;it<nt;it++){
				p1[isrc]+= -pow(dx,2.0)*ricker[it]*velextend[isrc];
		//printf("p2 no modelo:%f\n",p2[isrc]);	
	// Laplacian and field update
				for(int ix=(nderiv2-1);ix<(nxx-nderiv2);ix++){
					for(int iz=(nderiv2-1);iz<(nzz-nderiv2);iz++){
						gamma=gamma_x[ix] + gamma_z[iz];
						mgamma=-(1-gamma)/(1+gamma);
						invpgamma= 0.5*(1-mgamma);
						imodel=i0_model + iz + nzz*ix;
						laplacian=2*deriv2[0]*p2[imodel];
						for(int iconv=1;iconv < nderiv2;iconv++){
							laplacian+=deriv2[iconv]*p2[imodel - iconv];
							laplacian+=deriv2[iconv]*p2[imodel + iconv];
							laplacian+=deriv2[iconv]*p2[imodel - nzz*iconv];
							laplacian+=deriv2[iconv]*p2[imodel + nzz*iconv];
						}
						p1[imodel]= mgamma*p1[imodel] + invpgamma*( 2.0*p2[imodel] +
						velextend[imodel]*laplacian);
					//	p1[imodel]= -p1[imodel] + ( 2.0*p2[imodel] +
					//			velextend[imodel]*laplacian);
					}	
				}
				for(int iswap=0;iswap<nxx*nzz;iswap++){
					dswap=p2[iswap];
					p2[iswap]=p1[iswap];
					p1[iswap]=dswap;
				}
			modulo=fmod(it,ndt);
			if((int) modulo==0){
				im=0;
				for(int imx=0;imx<nx;imx++){
					for(int imz=0;imz<nz;imz++){
						//printf("im:%d\n",im);
						//printf("index p1: %d\n",i0_model + imz + imx*nzz);
						//printf("index p2: %d\n",im + iframe*nx*nz);
						//printf("value p1:%f\n",p1[im + iframe*nx*nz]);
						//printf("value pfield:%f\n",pfield[im + iframe*nx*nz]);
						pfield[im]=p1[i0_model + imz + imx*nzz];
						im+=1;
					}
				}
				fwrite(pfield,sizeof(float),nx*nz,wavemovie);
			}
	}
}
	printf("Tentando escrever");
//	int ntref;
//	float* paux;
//	FILE* waveim;
//	waveim=fopen("~/ProjetoFWI/models/waveimage.bin","wb");
//	paux=malloc(sizeof(float)*nx*nz);
//	/for(int ii=0;ii<nx*nz;ii++){
//		ntref=(int) 0.3/dtrec;
//		paux[ii]=pfield[ii + nx*nz*ntref];
//	}
//	fwrite(paux,sizeof(float),nx*nz,waveim);
	return 0;
}

