#include <stdio.h>      // Input and output
#include <stdlib.h>     // Pointers
#include <string.h>
#include <math.h>       // Math functions
#include <lapack.h>     // Linear algebra
#include "modelling_utils.h"
#define  PI 3.14159265358979323846264338328
#define dtrec 0.004
#define dtrtm 0.004
#define sec2ms 1.0e6
#define cfl 0.25
#define nderiv2 3
#define NB 200
void velextension(float* vex,float* v,int nx,int nz,struct b border,int nxx,int nzz){
	int i0;
	i0=(border.ixb+1)*nzz + border.izb+1;
	
	for(int ii=0;ii<nx;ii++){
            for(int jj=0;jj<nz;jj++){
                vex[(int) (i0 + ii*nzz + jj)]=v[ii*nz + jj];
            }
    }
    // Border at the left side
        for(int ii=0;ii<(border.ixb + 1);ii++){
            for(int jj=0;jj<nz;jj++){
					vex[border.izb +  1 + ii*nzz + jj]=v[jj];
             }
        }
        
    // Border at the right side
        for(int ii=0;ii<(border.ixb+1);ii++){
          for(int jj=0;jj<nz;jj++){
            vex[ (nzz*( border.ixb+1 +nx) +  border.izb + 1) + ii*nzz + jj]=v[jj];
          }
        }
   
	// Border at the top side
		for (int ii=0;ii<nxx;ii++){
			for(int jj=0;jj<(border.izb+1);jj++){
				vex[jj + ii*nzz] = vex[(border.izb+1) + ii*nzz];
			}
		}	
	// Border at the bottom side
		for (int ii=0;ii<nxx;ii++){
			for(int jj=border.ize;jj<nzz;jj++){
				vex[jj + ii*nzz] = vex[(border.ize-1) +  ii*nzz];
			}
		}	
}

double fricker(float t,float freq){
		double beta=PI*freq*t;
		double ricker;	
		ricker=(-sqrt(2.0)*beta + 1.0)*(1.0 + sqrt(2.0)*beta)*exp(-(beta*beta));
		return ricker;
}


int acoustic_modelling2D(char* velhdr,char* veldir,char* aquisitionhdr,char* secdir){
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
	float* dswap;
	int i0_model,imodel;
	float laplacian;


// Aquisition parameters
	int isrc,isx,isz,idxs,nshots;
	int it0,itrec;
	float xsmax,xsmin,dxs,zs,hmin,hmax,dh,zg ;
	float freq,ttotal;
	float* ricker;
	float* seismogram;
	
// Stability parameters
	int ns,nt,ndt;
	float maxv,minv;
	float dxmax,dt;
//	float deriv2[nderiv2]={-3.02359410, 1.75000000,-0.291666667, 0.0648148148,-0.0132575758,0.0021212122,-0.000226625227,0.0000118928690};
	float deriv2[nderiv2]={-5.0/2,4.0/3,-1.0/12};
	float deriv2_sum=0;

//Leitura do modelo de velocidade: hdr e bin
	int nparam,nvparam;
	FILE* VELHDR,*VELDIR,*PARAM,*SECTION;
    nparam=10;
	nvparam=6;
	VELHDR=fopen(velhdr,"r");  
    VELDIR=fopen(veldir,"rb");
	PARAM=fopen(aquisitionhdr,"r");

// Aquisiton parameters reading
	if (fscanf(PARAM,"%f %f %f %f %f %f %f %f %f %f",&xsmin,&xsmax,&dxs,&zs,&hmin,&hmax,&dh,&zg,&freq,&ttotal) < nparam){
			printf("Too few acquistion parameters!");
			return 0;
	}

// Velocity header reading
	if(fscanf(VELHDR,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx)< nvparam){
		printf("%d\n",nz);
		printf("%d\n",nx);
		printf("%f\n",z0);
		printf("%f\n",x0);
		printf("%f\n",dz);
		printf("%f\n",dx);
		printf("Missing velocity model header parameters!");
		return 0;
	}


	vel=malloc(sizeof(float)*nx*nz);
	if(vel==NULL){
		printf("Fail to allocate memory");
		exit(EXIT_FAILURE);
	}
	if(!fread(vel,sizeof(float),nx*nz,VELDIR)){
		printf("Error reading file!!\n");
		exit(EXIT_FAILURE);
	}
	fclose(VELDIR);
	fclose(VELHDR);
	fclose(PARAM);

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
	ricker=malloc(sizeof(double)*nt);

	
// Half duration of ricker number of indexes it0	
	it0 = floor(sqrt(6.0/PI)/(freq*dt)); 
	for(int it=0;it<nt;it++){
		ricker[it]=fricker((it-it0)*dt,freq);
	}
	if(dx > dxmax){
		printf("Frequency is too high");
		return 0;
	}

// Extending velocity model *******************************************
    border.nb=NB;
	nxx=2*border.nb + nx;
    nzz=2*border.nb + nz;
    border.ixb=border.nb -1;
	border.izb=border.ixb;
	border.ize=nzz - border.nb;
   	
	velextend=malloc(sizeof(float)*nxx*nzz);
	velextension(velextend,vel,nx,nz,border,nxx,nzz);
	
   
	printf("smin=%f\n",xsmin);
	printf("smax=%f\n",xsmax);
	printf("dxs=%f\n",dxs);
	printf("zs=%f\n",zs);
	printf("hmin=%f\n",hmin);
	printf("hmax=%f\n",hmax);
	printf("dh=%f\n",dh);
	printf("zg=%f\n",zg);
	printf("freq=%f\n",freq);
	printf("time=%f\n",ttotal);
	printf("Dt for stability: %.8f\n",dt);
	printf("Dx for stability:%f\n",dxmax);
	printf("ndt is:%i \n",ndt);
	printf("Maximum value of velocity: %f \n",maxv);
	printf("Minimum value of velocity: %f \n",minv);
	printf("Number of border points:%d\n",border.nb);



// Scaling  extended velocity for computational resource economy
	for(int ii=0;ii<nxx;ii++){
		for(int jj=0;jj<nzz;jj++){
			velextend[jj + ii*nzz]= pow(velextend[jj + ii*nzz],2.0)*pow(dt/dx,2.0);
		}
	}

// Wavefield variables 
	p1=calloc(nxx*nzz,sizeof(double));
	p2=calloc(nxx*nzz,sizeof(double));
    i0_model=nzz*(border.ixb+1) + (border.izb + 1); 
	printf("nt:%d\n",nt);

// Index of source in the model
	isx=(int) (xsmin/dx);  
	isz=(int) (zs/dz);	   
	isrc= i0_model + isz +  nzz*isx;
	idxs=(int)(dxs/dx);
	nshots=(int) ((xsmax - xsmin)/dxs);
	nshots+=1;

// Receivers parameters
	int ng,idxh;
	float modulo;
	idxh=(int) (dh/dx);
	ng=(int) ((hmax-hmin)/dh);
	ng+=1;
	float section[ns*ng];
	int igeo[ng];

// Absorption border parameters
	int ix,iz;
	float gamma_x[nxx];
	float gamma_z[nzz];
	float gamma, mgamma, invpgamma, beta;

	beta=PI*freq*dt;
	for(iz = 0; iz < nzz; ++iz) {
		gamma_z[iz]=0.0;
	}
	for (ix = 0; ix < nxx; ++ix){
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

	printf("nshots:%d\n",nshots);
	SECTION=fopen(secdir,"wb");
// Forward modelling --------------------------------------------------------- 
	for(int is=0;is<nshots;is++){
			printf("is:%d-%d\n",(is+1),nshots);
			itrec=0;
			for(int iswap=0;iswap<nxx*nzz;iswap++){
						p1[iswap]=0.0;
						p2[iswap]=0.0;
			}
			isrc=i0_model + isz +  nzz*isx + nzz*idxs*is;
			for(int ig=0;ig<ng;ig++){
				igeo[ig]=i0_model + (zg/dz) + nzz*(isx + idxs*is) + nzz*(hmin/dx)
						 + ig*idxh*nzz;
			}
		// Forward modelling --------------------------------------------------------- 
			for(int it=0;it<nt;++it){
					p1[isrc]+= -pow(dx,2.0)*ricker[it]*velextend[isrc];
		// Laplacian and field update
					for(ix=(nderiv2-1);ix<(nxx-nderiv2);ix++){
						for(iz=(nderiv2-1);iz<(nzz-nderiv2);iz++){
							gamma=gamma_x[ix] + gamma_z[iz];
							mgamma=-(1-gamma)/(1+gamma);
							invpgamma=0.5*(1-mgamma);
							imodel= iz + nzz*ix;
							laplacian=2*deriv2[0]*p2[imodel];
							for(int iconv=1;iconv < nderiv2;iconv++){
								laplacian+=deriv2[iconv]*p2[imodel - iconv];
								laplacian+=deriv2[iconv]*p2[imodel + iconv];
								laplacian+=deriv2[iconv]*p2[imodel - nzz*iconv];
								laplacian+=deriv2[iconv]*p2[imodel + nzz*iconv];
							}
							p1[imodel]= mgamma*p1[imodel] + invpgamma*(2.0*p2[imodel] +
							velextend[imodel]*laplacian);
						}	
					}
					
				dswap=p2;
				p2=p1;
				p1=dswap;
				modulo=fmod(it,ndt);
				if((int) modulo==0){
					// Section writing
					for(int ii=0;ii<ng;ii++){
						section[itrec + ns*ii]=p1[igeo[ii]];
					}
					itrec++;
				}
		}
		  fwrite(section,sizeof(float),ns*ng,SECTION);
	}
	return ns;
}







