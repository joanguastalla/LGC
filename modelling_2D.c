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
			/*
			// One-way border condition
	
				for (iz = 0;iz<nz;iz++) {
					for (ix = 0;ix<(nderiv2-1);ix++) {
						// Left x-border
						imodel=((nderiv2-2)-ix)*nzz + (border.nb + iz);
						p1[imodel]=(1-sqrt(velextend[imodel]))*p2[imodel] + 
	   					sqrt(velextend[imodel])*p2[imodel + nzz];
						// Right x-border
						imodel=((nxx-1)-ix)*nzz + (border.nb + iz);
						p1[imodel]=(1-sqrt(velextend[imodel]))*p2[imodel] +
	   					sqrt(velextend[imodel])*p2[imodel - nzz];
					}
				}
				for(ix=0;ix<nx;ix++){
					for(iz=0;iz<(nderiv2-1);iz++){
						// Top z-border
						imodel= (border.nb + ix)*nzz + iz;
						p1[imodel]=(1-sqrt(velextend[imodel]))*p2[imodel] +
	   					sqrt(velextend[imodel])*p2[imodel + 1];
						// Bottom z-border
						//imodel=(border.nb + ix)*nzz + ( (nzz-1) - (nderiv2-2) + iz);
						imodel=(border.nb + ix)*nzz + ( (nzz-1) - iz);
						p1[imodel]=(1-sqrt(velextend[imodel]))*p2[imodel] + 
	   					 sqrt(velextend[imodel])*p2[imodel-1];
					}
				}
			*/

			//modulo=fmod(it,ndt_rtm);
			/*
			precbackward[it]=p1[isrc];

			
			if((int) modulo==0){
						precbackward[itrec]=p1[isrc];
				itrec+=1;
			}
			*/
			 //pfield,p1(backward)  --> p
			/*
			 //Wavemovie writing
			if((int) modulo==0){
				im=0;
				for(int imx=0;imx<nx;imx++){
					for(int imz=0;imz<nz;imz++){
						p1aux[im]=p1[i0_model + imz + imx*nzz];
						im+=1;
					}
				}
				fwrite(p1aux,sizeof(float),nx*nz,wavemovie);
			}
			*/


	/*
	section=fopen(sectiondir,"wb");
	fwrite(precbackward,sizeof(float),nt*nshots,section);
	fclose(section);
	*/
/*
if(fwrite(prec,sizeof(float),ng*ns,section) < ng*ns){
	printf("Erro");
	exit(EXIT_FAILURE);
}
*/
//if(fwrite(pfield,sizeof(float),nx*nz*nt_rtm,wavemovie) < nx*nz*nt_rtm)
//	exit(EXIT_FAILURE);
int main(){

// Declaration of variables
    

// wavefield 
	float* p1,*p2;
	float* ref;
	float* pfield, *image,*p1aux;
	float* pforward;
	float* taper;
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
	float laplacian;
	float*	dswap;
// Aquisition parameters
	int it0,ng,itrec;
	int isrc,isx,isz,idxs,nshots;
	int* igeo;
	float* prec;
	float xsmax,xsmin,dxs,zs,hmin,hmax,dxh,zg;
	float freq,ttotal;
	double* ricker;
	float* sismogram;
	
// Stability parameters
	int ns,nt,nt_rtm,ndt,ndt_rtm,it_rtm=0;
	float maxv,minv;

	float dxmax,dt;
	//float deriv2[nderiv2]={-3.02359410, 1.75000000,-0.291666667, 0.0648148148,-0.0132575758,0.0021212122,-0.000226625227,0.0000118928690};
	float deriv2[nderiv2]={-5.0/2.0,4.0/3.0,-1.0/12.0};
	float deriv2_sum=0;


// Directories for reading or writing files  
    char* param_dir;
	char* velbin;
    char* velhdr;
	char* rickerdir;
    char* velextend_path;
	char* moviedir;
	char* sectiondir;
	char* imagedir;
//	param_dir="/export/home/joan/ProjetoFWI/param.txt" ;
	param_dir="/export/home/joan/ProjetoFWI/models/geometry.txt" ;
 //   velbin = "/export/home/joan/ProjetoFWI/models/refplane.bin";
//	velhdr = "/export/home/joan/ProjetoFWI/models/Vel.txt";
	velbin="/export/home/joan/ProjetoFWI/models/vel.bin";
	velhdr="/export/home/joan/ProjetoFWI/models/vel.hdr";
	velextend_path="./velextend.bin";
	rickerdir="/export/home/joan/ProjetoFWI/models/ricker.bin";
	moviedir="/scratch/joan/wavemovie.bin";
	sectiondir="/scratch/joan/section.bin";
	imagedir="/scratch/joan/image.bin";
// File pointers to open files
	FILE* velfile;
    FILE* velfile_hdr;
    FILE* velwrite;
    FILE* param;
	FILE* ricker_bin;	
    FILE* wavemovie;
	FILE* section;
	FILE* IMAGE;
	FILE* SISMOGRAM;
	// Alocação de matrizes e vetores
    
//Leitura do modelo de velocidade: hdr e bin
    velfile_hdr=fopen(velhdr,"r");    
    velfile=fopen(velbin,"rb");
	param=fopen(param_dir,"r");
	fscanf(param,"%f %f %f %f %f %f %f %f %f %f",&xsmin,&xsmax,&dxs,&zs,&hmin,&hmax,&dxh,&zg,&freq,&ttotal);
	fscanf(velfile_hdr,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx);
	ng=(int) (hmax - hmin)/dxh;
	ng+=1;
	igeo=malloc(sizeof(int)*ng);
	prec=calloc(ng,sizeof(float));
	vel=malloc(sizeof(float)*nx*nz);
	fread(vel,sizeof(float),nx*nz,velfile);
    fclose(velfile);
    fclose(velfile_hdr);
	fclose(param);

// Extending velocity model *******************************************
    border.nb=NB;
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
	printf("ns:%d\n",ns);
	dxmax=minv/(6.0*freq);
//	dt=dz/(4.0*maxv);
	for(int ii=0;ii<nderiv2;ii++){
		deriv2_sum+=abs(deriv2[ii]);
	}
	dt=cfl*sqrt(2.0/deriv2_sum)*dx/maxv;
	ndt=(int)(dtrec/dt);
	ndt+=1;
	ndt_rtm=(int)(dtrtm/dt);
	ndt_rtm+=1;
// Chosing dt for modelling to be a a divisor of dtrec
	if(ndt > 1){
		dt=(dtrec/ndt);
	}else{
		dt=dtrec;
	}
	nt=(int)(ttotal/dt);
	nt+=1;
	printf("nt:%d\n",nt);
	nt_rtm=(int)(ttotal/dtrtm);
	nt_rtm+=1;
	pforward=calloc(ng*nt,sizeof(float));
	ricker=malloc(sizeof(double)*nt);
	
	
	// Half duration of ricker number of indexes it0	
	it0 =4*floor(sqrt(6.0/PI)/(freq*dt)); 
	ricker_bin=fopen(rickerdir,"wb");
	for(int it=0;it<nt;it++){
		ricker[it]= fricker((it-it0)*dt,freq);
	}
	fwrite((float*) ricker,sizeof(float),nt,ricker_bin);
	fclose(ricker_bin);

    printf("xsmax:%.1f\n",xsmax);
    printf("xsmin:%.1f\n",xsmin);
    printf("dxs:%.1f\n",dxs);
    printf("zs:%.1f\n",zs);
    printf("hmin:%.1f\n",hmin);
    printf("hmax:%.1f\n",hmax);
    printf("dxg:%.1f\n",dxh);
    printf("zg:%.1f\n",zg);
    printf("freq:%.1f\n",freq);
    printf("ttotal:%.1f\n",ttotal);
	printf("dt|stable: %.8f\n",dt);
	printf("dx|stable:%.1f\n",dxmax);
	printf("nx:%d\n",nx);
	printf("nz:%d\n",nz);
	printf("dx:%.1f\n",dx);
	printf("dz:%.1f\n",dz);
	printf("ndt:%d \n",ndt);
	printf("max vel: %.1f \n",maxv);
	printf("min vel: %.1f \n",minv);
	printf("nzz:%d\n",nzz);
	printf("nxx:%d\n",nxx);
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

	
// Wavefield temporal evolution common-shot situation
	p1=calloc(nxx*nzz,sizeof(float));
	p2=calloc(nxx*nzz,sizeof(float));
    i0_model=nzz*(border.ixb+1) + (border.izb + 1); 
	printf("i0model:%d\n",i0_model);
// Index of source in the model
	
	isx=(int) (xsmin/dx);   // starting at 0 sample
	isz=(int) (zs/dz);	   // starting at 0 sample
	isrc= i0_model + isz +  nzz*isx;
	idxs=(int)(dxs/dx);
	nshots=(int) ((xsmax - xsmin)/dxs);
	nshots+=1;
// Wave movie variables 
ntrec=(int) (nt/ndt);
pfield=calloc(nx*nz*nt_rtm,sizeof(float));
p1aux=calloc(nx*nz,sizeof(float));
image=calloc(nx*nz,sizeof(float));
volatile float modulo;
// Wave propagation in borders
int ix,iz,idxh;
idxh=(int) (dxh/dx);
float* gamma_x;
float* gamma_z;
float gamma, mgamma, invpgamma, beta;
gamma_x=malloc(sizeof(float)*nxx);
gamma_z=malloc(sizeof(float)*nzz);

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


sismogram=calloc(sizeof(float),ns*ng);
memset(image,0,sizeof(float)*nx*nz);
wavemovie=fopen(moviedir,"wb");
SISMOGRAM=fopen("/scratch/joan/cshot.bin","wb");
// Forward modelling --------------------------------------------------------- 
for(int is=0;is<nshots;is++){
		printf("is:%d-%d\n",is,nshots);
		itrec=it_rtm=0;
		memset(pfield,0,nx*nz*nt_rtm*sizeof(float));
		for(int iswap=0;iswap<nxx*nzz;iswap++){
					p1[iswap]=0.0;
					p2[iswap]=0.0;
		}
		isrc=i0_model + isz +  nzz*isx + nzz*idxs*is;
		printf("isrc:%d\n",isrc);
		iframe=0;
		for(int ig=0;ig<ng;ig++){
			igeo[ig]=i0_model + (zg/dz) + nzz*(isx + idxs*is) + nzz*(hmin/dx)
					 + ig*idxh*nzz;
		}
	// Forward modelling --------------------------------------------------------- 
		for(int it=0;it<nt;++it){
				p1[isrc]+= -pow(dx,2.0)*ricker[it]*velextend[isrc];
	// Laplacian and field update
				for(int ix=(nderiv2-1);ix<(nxx-nderiv2);ix++){
					for(int iz=(nderiv2-1);iz<(nzz-nderiv2);iz++){
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
			
			for(int ig=0;ig < ng;ig++){
				pforward[it + nt*ig]=p1[igeo[ig]];
			}
			modulo=fmod(it,ndt_rtm);
			if((int) modulo==0){
				// Wavemovie writing
				im=0;
				for (int ig = 0;ig < ng;++ig){
					sismogram[it_rtm + ig*ns]=p1[igeo[ig]];
				}
				for(int imx=0;imx<nx;imx++){
					for(int imz=0;imz<nz;imz++){
						p1aux[im]=p1[i0_model + imz + imx*nzz];
						im+=1;
					}
				}
				fwrite(p1aux,sizeof(float),nx*nz,wavemovie);
				memcpy((pfield + it_rtm*nx*nz),p1aux,sizeof(float)*nx*nz);
				it_rtm++;
			}
	}
	fwrite(sismogram,sizeof(float),ns*ng,SISMOGRAM);
	memset(p1,0,sizeof(float)*nxx*nzz);
	memset(p2,0,sizeof(float)*nxx*nzz);
	it_rtm=0;
	itrec=0;
	// Backward modelling
	for(int it=0;it < nt;it++){
		for(int is=0;is<ng;is++){ 
			iframe=0;
			p1[igeo[is]]+= -pow(dx,2.0)*pforward[(nt-1) - it + is*nt]*velextend[igeo[is]];
		}
		// Laplacian and field update
		for(int ix=(nderiv2-1);ix<(nxx-nderiv2);ix++){
			for(int iz=(nderiv2-1);iz<(nzz-nderiv2);iz++){
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
		modulo=fmod(it,ndt_rtm);
		dswap=p2;
		p2=p1;
		p1=dswap;
		if((int) modulo==0){
			for(int imx=0;imx<nx;++imx){
				for(int imz=0;imz<nz;++imz){
					imodel=imz + imx*nz;
					image[imodel]+=pfield[imodel + ((nt_rtm-1)-it_rtm)*nx*nz]*
					p1[i0_model + imz + imx*nzz];
				}
			}
			it_rtm++;
		}
			
	}
}




ref=malloc(sizeof(float)*nx*nz);
// Laplacian filter on image
for(int imx=(nderiv2-1);imx < (nx-nderiv2); ++imx){
	for(int imz=(nderiv2-1);imz<(nz-nderiv2); ++imz){
				imodel=imz + nz*imx;
				laplacian=2*deriv2[0]*image[imodel];
				for(int iconv=1;iconv < nderiv2;iconv++){
					laplacian+=deriv2[iconv]*image[imodel - iconv];
					laplacian+=deriv2[iconv]*image[imodel + iconv];
					laplacian+=deriv2[iconv]*image[imodel - nz*iconv];
					laplacian+=deriv2[iconv]*image[imodel + nz*iconv];
				}
				ref[imodel]=laplacian;
	}
}
// Tappering
taper=calloc(nz,sizeof(float));
int taper_end,taper_begin,ntaper;
FILE* TAPER;
char* taper_dir="./taper.bin";
taper_begin=21;
taper_end=51;
ntaper=taper_end-taper_begin;
for(int ii=taper_end;ii<nz;ii++){
	taper[ii]=1.0;
}

for(int ii=0;ii<ntaper;ii++){
	taper[taper_begin + ii]=exp((ii-ntaper)/dz);
}
for(int ii=0;ii<nz;ii++){
	printf("taper[%d]=%f\n",ii,taper[ii]);
}


for (int imx=0;imx<nx;++imx) {
	for (int imz=0;imz< nz; ++imz){
		ref[imx*nz + imz]=ref[imx*nz + imz]*taper[imz];
	}	
}

TAPER=fopen(taper_dir,"wb");
fwrite(taper,sizeof(float),nz,TAPER);


free(p1);
free(p2);
fclose(wavemovie);
IMAGE=fopen(imagedir,"wb");
fwrite(ref,sizeof(float),nx*nz,IMAGE);
fclose(IMAGE);
return 0;
}

