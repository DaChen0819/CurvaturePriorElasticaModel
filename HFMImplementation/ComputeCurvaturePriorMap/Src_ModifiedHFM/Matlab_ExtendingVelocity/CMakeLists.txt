cmake_minimum_required(VERSION 3.6)
set (CMAKE_CXX_STANDARD 17)

project(HamiltonExtendingVelocity)

# Matlab setting
set(MATLAB_APP "/Applications/MATLAB_R2023a.app")
set(MATLAB_LIB_SUBDIR "/bin/maci64")
set(MATLAB_INCLUDE_SUBDIR "/extern/include")

set(CMAKE_MACOSX_RPATH 0)

set(JMM_CPPLibs_dir "${CMAKE_CURRENT_SOURCE_DIR}/../JMM_CPPLibs" CACHE FILEPATH "JMM_CPPLibs directory")
set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")
add_subdirectory("../JMM_CPPLibs" "${DummyBinDir}/JMM_CPPLibs")
add_subdirectory("../GeoMetrics" "${DummyBinDir}/GeoMetrics")
add_subdirectory("FastMarchingBase" "${DummyBinDir}/FastMarchingBase")

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
include_directories("${MATLAB_APP}${MATLAB_INCLUDE_SUBDIR}")
include_directories("../MatlabInterface")

#cmake_policy(SET CMP0057 NEW)

set(MATLAB_LIB_DIR "${MATLAB_APP}${MATLAB_LIB_SUBDIR}")
set(MATLAB_LIBRARIES "${MATLAB_LIB_DIR}/libmat.dylib" "${MATLAB_LIB_DIR}/libmex.dylib" "${MATLAB_LIB_DIR}/libmx.dylib")


set(Project_Headers 
 	${Base_Headers}
 	${Base_Implementation_Headers}
 	${Specializations_Headers}
 	${Specializations_Implementation_Headers}
 	${Experimental_Headers}
 	${Experimental_Implementation_Headers}
 	${ExtraAlgorithms_Headers}
 	${ExtraAlgorithms_Implementation_Headers}
 	${Root_Headers}
 	${DataStructures_Headers}
 	${LinearAlgebra_Headers}
 	${Output_Headers}
 	${Macros_Headers}
)


 # ---- Headers IDE layout -----
source_group("DataStructures" FILES ${DataStructures_Headers})
source_group("Macros" FILES ${Macros_Headers})
source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
source_group("Output" FILES ${Output_Headers})
source_group("Base" FILES ${Base_Headers})
source_group("Base\\Implementation" FILES ${Base_Implementation_Headers})
source_group("Specializations" FILES ${Specializations_Headers})
source_group("Specializations\\Implementation" FILES ${Specializations_Implementation_Headers})
source_group("Experimental" FILES ${Experimental_Headers})
source_group("Experimental\\Implementation" FILES ${Experimental_Implementation_Headers})
source_group("ExtraAlgorithms" FILES ${ExtraAlgorithms_Headers})
source_group("ExtraAlgorithms\\Implementation" FILES ${ExtraAlgorithms_Implementation_Headers})
source_group("RootInterface" FILES  ${Root_Headers})

add_library(HamiltonExtendingVelocity SHARED HamiltonExtendingVelocity.cxx ${Root_Headers} ${Base_Headers} ${Base_Implementation_Headers} ${ExtraAlgorithms_Headers} ${ExtraAlgorithms_Implementation_Headers} ${Specializations_Headers} ${Specializations_Implementation_Headers} ${Experimental_Headers} ${Experimental_Implementation_Headers})

target_link_libraries(HamiltonExtendingVelocity ${MATLAB_LIBRARIES})









