add_subdirectory(Specializations)
set(Specializations_Headers ${Specializations_Headers} PARENT_SCOPE)
set(Specializations_Implementation_Headers ${Specializations_Implementation_Headers} PARENT_SCOPE)


add_subdirectory(Experimental)
set(Experimental_Headers ${Experimental_Headers} PARENT_SCOPE)
set(Experimental_Implementation_Headers ${Experimental_Implementation_Headers} PARENT_SCOPE)

set(PROJECT_DIR ${PROJECT_DIR} "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
