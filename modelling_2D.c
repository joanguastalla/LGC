#include <stdio.h>      // Input and output
#include <stdlib.h>     // Pointers
#include <string.h>
#include <math.h>       // Math functions
#include <lapack.h>     // Linear algebra
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
	int nxx,nzz,nb;
    int i0;
    int ixb,izb,ixe,ize;
	float* velextend;
	
	// Parametros de aquisição
	float xfmax,xfmin,dxf,zf,xg_offset_min,xg_offset_max,dxg,zg ;
	float freq,ttotal;
	
	// Parametros de estabilidade
	int ns,ndt;
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
        nb=50;
        ixb=nb-1;
        izb=ixb;
		nxx=2*nb + nx;
        nzz=2*nb + nz;
        ize=nzz - nb;
		i0=(ixb+1)*nzz + izb+1;
        velextend=malloc(sizeof(float)*nxx*nzz);
   
   // First loop for velocity extending
    
   
    // Apparently a allocated pointer intializes to zero in all positions, 
	// when no values are passed in

    for(int ii=0;ii<nx;ii++){
            for(int jj=0;jj<nz;jj++){
                velextend[(int) (i0 + ii*nzz + jj)]=vel[ii*nz + jj];
            }
    }
    // Border at the left side
        for(int ii=0;ii<(ixb + 1);ii++){
            for(int jj=0;jj<nz;jj++){
                velextend[izb +  1 + ii*nzz + jj]=vel[jj];
             }
        }
        
    // Border at the right side
        for(int ii=0;ii<(ixb+1);ii++){
            for(int jj=0;jj<nz;jj++){
                velextend[ (nzz*( ixb+1 +nx) +  izb + 1) + ii*nzz + jj]=vel[jj];
             }
        }
   
	// Border at the top side
		for (int ii=0;ii<nxx;ii++){
			for(int jj=0;jj<(izb+1);jj++){
				velextend[jj + ii*nzz] = velextend[(izb+1) + ii*nzz];
			}
		}	
	// Border at the bottom side
		for (int ii=0;ii<nxx;ii++){
			for(int jj=ize;jj<nzz;jj++){
				velextend[jj + ii*nzz] = velextend[(ize-1) +  ii*nzz];
			}
		}	
   
	
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

	ns=(int) (ttotal/dtrec) ;
	ns+=1;
	dxmax=minv/(6.0*freq);
	for(int ii=0;ii<nderiv2;ii++){
		deriv2_sum+=abs(deriv2[ii]);
	}
	dt=cfl*sqrt(2.0/deriv2_sum)*dx/maxv;
	ndt=(int) (dtrec/dt);
	ndt+=1;
	printf("Dt for stability: %f\n",dt);
	printf("Maximum value of velocity: %f \n",maxv);
	printf("Minimum value of velocity: %f \n",minv);
	return 0;
	


}

