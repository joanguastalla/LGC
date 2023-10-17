#include <math.h>
#define  PI 3.14159265358979323846264338328
struct b{
	int ixb;
	int izb;
	int ize;
	int nb;
};
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

float fricker(float t,float freq){
		float beta=pow(PI*freq*t,2);
		float ricker;	
		ricker=(1.0 - 2*beta)*exp(-beta);
		return ricker;
}
