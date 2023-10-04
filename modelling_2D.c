#include <stdio.h>      // Input and output
#include <stdlib.h>     // Pointers
#include <string.h>
#include <math.h>       // Math functions
#include <lapack.h>     // Linear algebra

int main(){
       
// Declaracao de Variaveis
    // Parametros do modelo original e extendido
    int nx,nz,nxx,nzz,nb;
    int i0;
    int ixb,izb,ixe,ize;
    float dx,dz,x0,z0;
    float* vel;
    float* velextend;
// Ler o modelo de velocidade
    FILE* velfile;
    FILE* velfile_hdr;
    FILE* velwrite;
    char* velbin;
    char* velhdr;
    char* velextend_path;
    // Alocação de matrizes e vetores
    
    //Leitura do modelo de velocidade: hdr e bin
    velhdr = "/export/home/joan/ProjetoFWI/models/Vel.txt";
    velbin = "/export/home/joan/ProjetoFWI/models/Vel.bin";
    velfile_hdr=fopen(velhdr,"r");    
    velfile=fopen(velbin,"rb");

    fscanf(velfile_hdr,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx);
    vel=malloc(sizeof(float)*nx*nz);
	fread(vel,sizeof(float),nx*nz,velfile);
    fclose(velfile);
    fclose(velfile_hdr);
    //printf("nz= %i nx= %i z0=%f x0=%f dz=%f dx=%f \n",nz,nx,z0,x0,dz,dx);
    //printf("v0=%f v1=%f ", vel[1], vel[500]);   
  
    // Extending velocity model
        nb=20;
        ixb=nb-1;
        izb=ixb;
        nxx=2*nb + nx;
        nzz=2*nb + nz;
        i0=(ixb+1)*nzz + izb+1;
        velextend=malloc(sizeof(float)*nxx*nzz);
   
   // First loop for velocity extending
    
    printf("i0: %i \n",i0);
   
    // Apparently a allocated pointer intializes to zero in all positions, when no values are passed in

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
    
    printf("first_extend: %f \n",velextend[i0+10]);
    velextend_path="./velextend.bin";
    velwrite=fopen(velextend_path,"wb");
    fwrite(velextend,4,nxx*nzz,velwrite);    
    
    return 0;
}

