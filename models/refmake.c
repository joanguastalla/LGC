#include <stdio.h>

void main(){
	int nx,nz,idep;
	float dx,dz,depth,v0,v1;
	FILE* VMODEL;
	char* vmodel="refplane.bin";
	nx=401;
	nz=301;
	dx=dz=4;
	depth=600;
	v0=2000;
	v1=2400;
	idep=(depth/dz) + 1;
	printf("idep:%d\n",idep);
	float V[nx][nz];
	for(int ii=0;ii<nx;++ii){
		for(int jj=0;jj<idep;jj++){
			V[ii][jj]=v0;
		}
	}
	for(int ii=0;ii<nx;++ii){
		for(int jj=idep;jj<nz;jj++){
			V[ii][jj]=v1;
		}
	}
	VMODEL=fopen(vmodel,"wb");
	fwrite(V,sizeof(float),nx*nz,VMODEL);

}
