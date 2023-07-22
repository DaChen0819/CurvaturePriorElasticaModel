#ifndef DispatchAndRun_h
#define DispatchAndRun_h
 
#include "JMM_CPPLibs/Macros/String.h"
#include "JMM_CPPLibs/Macros/PPCat.h"

#include "Specializations/Curvature2.h"
#include "Specializations/Curvature3.h"
#include "Experimental/PrescribedCurvature2.h"
#include "Experimental/RegionalCurvature2.h"

// ------- Custom invocation, with multiple models.  ---------
#define HFMSpecializationMacro(modelName) \
{ \
using StencilDataType = Stencil ## modelName ;\
using HFMI = StencilDataType::HFMI; \
if(model== #modelName){ \
    io.currentSetter=IO::SetterTag::Compute;\
    StencilDataType stencil; \
    HFMI(io, stencil).Run();\
    io.currentSetter=IO::SetterTag::User; return;} \
}

void Run(IO & io){
    const std::string rawModel = io.GetString("model");
    if(rawModel == "Elastica2"){
        const std::string model="Elastica2<5>";
        HFMSpecializationMacro(Elastica2<5>);
    }
    else if(rawModel == "Elastica2_9"){
        const std::string model="Elastica2<9>";
        HFMSpecializationMacro(Elastica2<9>);
    }
    else if(rawModel == "ReedsShepp2"){
        const std::string model="ReedsShepp2";
        HFMSpecializationMacro(ReedsShepp2);
    }
    else if(rawModel == "ReedsSheppExt2"){
        const std::string model="ReedsSheppExt2";
        HFMSpecializationMacro(ReedsSheppExt2);
    }
    else if(rawModel == "ReedsSheppForward2"){
        const std::string model="ReedsSheppForward2";
        HFMSpecializationMacro(ReedsSheppForward2);
    }
    else if(rawModel == "Dubins2"){
        const std::string model="Dubins2";
        HFMSpecializationMacro(Dubins2);
    }
    else{
        ExceptionMacro("Unrecognized model : " << rawModel);
    }
    
}


#endif 
