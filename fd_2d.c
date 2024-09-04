#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "modelling_utils.h"
#define  PI 3.14159265358979323846264338328
#define nderiv2 3
#define dtrec 4e-3
#define sec2ms 1.0e6
#define cfl 0.25
#define NB 100
void print_help() {
	printf("FD2D - Finite differences (4th order) acoustic wave modelling.\n\n");
    printf("Usage: fd2d < vhdr vbin geom stdout  [options]\n");
	printf("\n");
    printf("Required parameters:\n\n");
	printf("vhdr - Velocity model header (nz nx x0 z0 dx dz)\n");
	printf("vbin - Binary velocity model\n");
	printf("geom - CS Aquisiton geometry (xsmin  xsmax  dxs  zs  h0  hf  dh  zg  freq  ttotal) \n");
	printf("stdout - Filepath for saving modelled data\n\n");
	printf("--------------------------------------------------\n");
	printf("nz -> number of samples in depth\n");
	printf("nx -> number of horizontal samples\n");
	printf("x0 -> First x coordinate in velocity model\n");
	printf("z0 -> First z coordinate in velocity model\n");
	printf("dx -> sampling in x-direction\n");
	printf("dx -> sampling in z-direction\n");
	printf("--------------------------------------------------\n");
	printf("--------------------------------------------------\n");
	printf("xsmin -> initial shotpoint\n");
	printf("xsmax -> final shotpoint\n");
	printf("dxs -> shots sampling\n");
	printf("zs -> source depth\n");
	printf("h0 -> initial offset\n");
	printf("hf -> final offset\n");
	printf("dh -> offset sampling\n");
	printf("zg -> receiver depth\n");
	printf("freq -> peak frequency of Ricker wavelet\n");
	printf("ttotal ->  total time of wave propagation\n");
	printf("--------------------------------------------------\n");
	printf("Options:\n");
    printf("  -h, --help       Display this help message and exit\n");
}
int main(int argc,char** argv){
	char* veldir="/export/home/joan/Seismic-Datasets/BP/vp.rsf@";
	char* velhdr="/export/home/joan/Seismic-Datasets/BP/vel_hdr.txt";
	char* paramdir="/export/home/joan/Seismic-Datasets/BP/geometry.txt";
	char* secdir="/scratch/joan/bp_cshot.bin";
	for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0){
            print_help();
            return 0;
        }
    }
	int	nvar_funcall=4;
	if((argc-1)!=nvar_funcall){
		printf("You must provide velocity,velocity header, aquisition header and section saving dir i				n that order!!\n");
		exit(EXIT_FAILURE);
	}
	acoustic_modelling2D(argv[1],argv[2],argv[3],argv[4]);
	return 0;
}




