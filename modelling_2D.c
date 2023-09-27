#include <stdio.h>      // Input and output
#include <stdlib.h>     // Pointers
#include <string.h>
#include <math.h>       // Math functions
//#include <glib2.h>
//#include <lapack.h>     // Linear algebra

int main(){
       
// Declaracao de Variaveis
    // Parametros do modelo original e extendido
    int nx,nz,nxx,nzz;
    int ixb,izb,ixe,ize;
    float dx,dz,x0,z0;
    float* vel;

// Ler o modelo de velocidade
    FILE* velfile;
    FILE* velfile_hdr;
    FILE* velwrite;
    char* velbin;
    char* velhdr;
    // Alocação de matrizes e vetores
    
    //Leitura do modelo de velocidade: hdr e bin
    velhdr = "/export/home/andres/Mestrado/FINAL/Dos_camadas/Caso1/velmodhom.hdr";
    velbin = "/export/home/andres/Mestrado/FINAL/Dos_camadas/Caso1/velhomo1.bin";
    velfile_hdr=fopen(velhdr,"r");    
    velfile=fopen(velbin,"rb");

    fscanf(velfile_hdr,"%i %i %f %f %f %f",&nz,&nx,&z0,&x0,&dz,&dx);
    vel=malloc(sizeof(float)*nx*nz);
	fread(vel,sizeof(float),nx*nz,velfile);
    fclose(velfile);
    fclose(velfile_hdr);
    printf("nz= %i nx= %i z0=%f x0=%f dz=%f dx=%f \n",nz,nx,z0,x0,dz,dx);
    printf("v0=%f v1=%f ", vel[1], vel[2]);
    
    return 0;
}

