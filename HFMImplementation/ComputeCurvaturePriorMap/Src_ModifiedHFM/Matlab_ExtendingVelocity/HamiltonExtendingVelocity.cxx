#include "JMM_CPPLibs/Output/MexIO.h"
typedef IO_<MexIO> IO;
typedef typename IO::Msg Msg;
typedef typename IO::WarnMsg WarnMsg;

#include "FastMarchingBase/DispatchAndRun.h"

void mexFunction(int nlhs, mxArray *plhs[],
                                     int nrhs, const mxArray *prhs[] ){
    if(nrhs!=1 || nlhs!=1){
        IO::WarnMsg() << "Exactly one input and one output are expected (in structures).\n";
        return;
    }
    
    
    try{
        IO io(prhs[0],plhs);
        io.arrayOrdering = IO::ArrayOrdering::YXZ_ColumnMajor;
        Run(io);
    }
    catch(const std::exception & e){
        IO::WarnMsg() << "Hamiltonian Fast Marching exception.\n " << e.what() << "\n";
    }
}

