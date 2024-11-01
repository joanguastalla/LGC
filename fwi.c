// Modelling data with smoothed model
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modelling_utils.h"

int main(){
		
	// AQUISITION PARAMETERS
	float xsmin,xsmax,dxs,zs,hmin,hmax,dh,zg,freq,ttotal,nparam=10;
	float x0,z0,dz,dx;
	int nz,nx;
	int ng,ns,nshots;
	FILE *OBS,*CALC,*GEOMETRY,*GRAD,*PARAM,*VELHDR,*VEL;
	float *obs,*calc,*residuo,*grad,*gradshot,*vel,*velaux;
	char* velhdr_obs="./models/vel.hdr";
	char* veldir_obs="./models/vel.bin";
	char* param_obs="./models/geometry.txt";
	char* secdir_obs="./models/cshot_obs.bin";
	char* secdir_observed="./models/cshot_obs.bin";
	// Optimization parameters
	int nepochs=50;
	float step=1.0;
	float *loss=calloc(sizeof(float),nepochs),loss_step,maxgrad;
	FILE *LOSS=fopen("loss.bin","wb");
	PARAM=fopen(param_obs,"r");
	VELHDR=fopen(velhdr_obs,"r");
	// Aquisiton parameters reading
	if (fscanf(PARAM,"%f %f %f %f %f %f %f %f %f %f",&xsmin,&xsmax,&dxs,&zs,&hmin,&hmax,&dh,&zg,&freq,&ttotal) < nparam){
			printf("Too few acquistion parameters!");
			exit(EXIT_FAILURE) ;
	}
	GEOMETRY=fopen(param_obs,"r");
	fscanf(GEOMETRY,"%f %f %f %f %f %f %f %f %f %f",&xsmin,&xsmax,&dxs,&zs,&hmin,&hmax,&dh,&zg,&freq,&ttotal);
	fscanf(VELHDR,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx);
	//  OBSERVED DATA ON EXACT VELOCITY MODEL
	fclose(VELHDR);
	fclose(GEOMETRY);
	ns=acoustic_modelling2D(velhdr_obs,veldir_obs,param_obs,secdir_obs);
	ng=(hmax-hmin)/dh;
	ng++;
	nshots=(xsmax-xsmin)/dxs;
	nshots++;
	//MODELLING ON INEXACT VELOCITY MODEL
	char* velhdr_calc="./models/vel.hdr";
	char* veldir_calc="./models/refplane_smooth.bin";
	char* velfwi_calc="./models/velfwi.bin";
	char* param_calc="./models/geometry.txt";
	//char* paramtemp="./models/geometry_header.txt";
	char* secdir_calc="./models/cshot_calc.bin";
	char cp[200];
	grad=calloc(sizeof(float),nx*nz);
	obs=calloc(sizeof(float),ns*ng);
	calc=malloc(sizeof(float)*ns*ng);
	residuo=malloc(sizeof(float)*ns*ng);
	VEL=fopen(velfwi_calc,"rb");
	vel=malloc(sizeof(float)*nx*nz);
	velaux=malloc(sizeof(float)*nx*nz);
	fread(vel,sizeof(float),nx*nz,VEL);
	fclose(VEL);
	// Command to copy initial velfwi velocity
	//sprintf(cp,"cp %s %s",veldir_calc,velfwi_calc);
	//printf("Result of command execution:%d\n",system(cp));
	printf("------------------- Start of RTM Gradient calculation ------------------------\n\n");
	for(int iepoch=0;iepoch<nepochs;iepoch++){
		maxgrad=0;
		OBS=fopen(secdir_obs,"rb");
		CALC=fopen(secdir_calc,"rb");
		for(int ifonte=0;ifonte<nshots;ifonte++){
			printf("------------------- RTM Gradient: shot %d\n",ifonte);
			acoustic_modelling2Dshot(velhdr_calc,velfwi_calc,param_calc,secdir_calc,ifonte);
			fread(obs,sizeof(float),ns*ng,OBS);
			fread(calc,sizeof(float),ns*ng,CALC);
			for(int ig = 0; ig < ng; ++ig){
				for(int ii = 0; ii < ns; ++ii){
					residuo[ii + ns*ig]=obs[ii + ns*ig] - calc[ii + ns*ig];
					loss[iepoch]+=residuo[ii + ns*ig]*residuo[ii + ns*ig];
				}
			}
			gradshot=rtm(velhdr_calc,velfwi_calc,param_calc,secdir_calc,residuo,ifonte);
			for(int ix = 0; ix < nx; ++ix){
				for(int iz = 0; iz < nz; ++iz){
					grad[ix*nz + iz]+=gradshot[ix*nz + iz];
				}
			}
			free(gradshot);
		}
		fclose(OBS);
		fclose(CALC);
		// Tapering source gradients
		for(int imz=0;imz<10;imz++){
			for(int imx=0;imx<nx;imx++){
				grad[imx*nz + imz]=0.0;
			}
		}
		for(int im=0;im<nx*nz;im++){
			if(grad[im]>maxgrad){
				maxgrad=grad[im];
			}
		}
		for(int im=0;im<nx*nz;im++){
			grad[im]/=maxgrad;
			grad[im]*=50;
		}
		loss[iepoch]/=2;
		// Optimzation Step
		GRAD=fopen("./grad.bin","wb");
		fwrite(grad,sizeof(float),nx*nz,GRAD);
		fclose(GRAD);
		printf("loss_epoch:%f\n",loss[iepoch]);
		// Loop for determining step size
		while(1){
			loss_step=0;
			for(int ii=0;ii<nx*nz;ii++){
					velaux[ii]=vel[ii] + grad[ii]*step;
			}
			VEL=fopen(velfwi_calc,"wb");
			fwrite(velaux,sizeof(float),nx*nz,VEL);
			fclose(VEL);
			OBS=fopen(secdir_obs,"rb");
			CALC=fopen(secdir_calc,"rb");
			printf("Half cutting the step length\n");
			for(int ifonte=0;ifonte<nshots;ifonte++){
				acoustic_modelling2Dshot(velhdr_calc,velfwi_calc,param_calc,secdir_calc,ifonte);
				fread(obs,sizeof(float),ns*ng,OBS);
				fread(calc,sizeof(float),ns*ng,CALC);
				for(int ig = 0; ig < ng; ++ig){
					for(int ii = 0; ii < ns; ++ii){
						residuo[ii + ns*ig]=obs[ii + ns*ig] - calc[ii + ns*ig];
						loss_step+=residuo[ii + ns*ig]*residuo[ii + ns*ig];
					}
				}
			}
			fclose(OBS);
			fclose(CALC);
			loss_step/=2;
			if(loss[iepoch]>loss_step){
					break;
			}
			step/=2;
			printf("Loss:%f\n",loss_step);
		}
		for(int iv=0;iv<nx*nz;iv++){
			vel[iv]=velaux[iv];
		}
	}
	fwrite(loss,sizeof(float),nepochs,LOSS);
	return 0;
}
