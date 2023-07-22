#include "mex.h"
#include <cmath>
#include <iostream>
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrive arguments */
	if( nrhs!=2 )
		mexErrMsgTxt("2 input arguments are required.");
	if( nlhs!=1 )
		mexErrMsgTxt("1 output arguments are required.");
    
    
	size_t nx = mxGetM(prhs[0]);
	size_t ny = mxGetN(prhs[0]);
    
	double* skeletonImage =   mxGetPr(prhs[0]);
    double* branchPtImage =   mxGetPr(prhs[1]);
    
	// first ouput : first component
	plhs[0] = mxCreateDoubleMatrix(nx,ny, mxREAL);
	double* output = mxGetPr(plhs[0]);
    
// some definitions.
#define ACCESS_ARRAY(A,i,j)   A[(i)+nx*(j)]
#define skeletonImage_(i,j)   ACCESS_ARRAY(skeletonImage,i,j)
#define branchPtImage_(i,j)   ACCESS_ARRAY(branchPtImage,i,j)
#define output_(i,j)          ACCESS_ARRAY(output,i,j)

    for(int i=0;i<nx;i++){
        for(int j=0; j<ny;j++){
            output_(i,j)=skeletonImage_(i,j);
        }
    }
    
    for(int i=0;i<nx;i++){
        for(int j=0; j<ny;j++){
            double isBranch=branchPtImage_(i,j);
            if(isBranch>0.5){
                for(int a = -1;a<=1;a++){
                    for(int b = -1;b<=1;b++){
                        if(i+a>=0 && i+a<nx && j+b<ny && j+b>=0){
                            int index=i+a+((int)nx)*(j+b);
                            output[index]=0.0;
                        }
                    }
                }
            }
        }
    }

	return;
}




