cmake_minimum_required(VERSION 2.6)

project(GMIN)

set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "dw34@cam.ac.uk")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)
enable_language(Fortran)

#SET (CMAKE_Fortran_FLAGS_RELEASE "-ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv")
SET (CMAKE_Fortran_FLAGS_RELEASE "-Mextend -O3 -Munroll -Mnoframe")

#include(CheckCXXCompilerFlag)
#find_package(LAPACK REQUIRED) 
file(GLOB DUMMY_AMBER amhdummy_cmake.f)
file(GLOB DUMMY_CHARMM chdummy.f)
file(GLOB DUMMY_AMBER9 amber9dummy.f90)

#file(GLOB GMIN_SOURCES1 porfuncs.f90 commons.f90 countatoms.f90 modamb.f90 modcharmm.f90 modmxatms.f90 modhess.f90 modamber9.f90 grouprotation.f90)
file(GLOB GMIN_SOURCES *.f *.f90 AMH/amhglobals.f AMH/amh_interfaces.f90 AMH/E_write.f90)
file(GLOB NOT_GMIN_SOURCES 
    molinfo.f90 tetherfuncall.f90 transforms.f90 olami.f rdpot.f wenzel.1D.f wenzel.2D.f compress2.f 
    mf.f otp.f bmin.f evstep.f lbfgs.f takestep2.f ptbs.f io_sanity.f90 fsc1TH.f90 dblpy.f90 funcs.f90
    dfuncs.f90 io.f90 mek-quake.f90 wham_djw.f90 dgetf2.f  dtrti2.f charmmBildc.f blas.f
    lapack.f dtrtri.f mind.f nonATLASblas.f dlaswp.f  symmetry3.f90 bipartite.f90    
    dlamch.f myblas.f mylapack.f rk45ad.f symmetry2.f90 tn.f amberdummy.f
    mc.f bspt.f display_version.f90 porfuncs.f90 main.f # remove these in calse already exist. will be added later
)
list(REMOVE_ITEM GMIN_SOURCES ${NOT_GMIN_SOURCES} ${DUMMY_AMBER} ${DUMMY_CHARMM} ${DUMMY_AMBER9})

################ begin More elegany solution needed! ####################
# these custom commands are to do the preprocessing of fotran files
# there is a nicer way, have to figure out how
add_custom_command (OUTPUT mc.f DEPENDS mc.F COMMAND cpp -traditional -P mc.F > mc.f)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES mc.f)

add_custom_command (OUTPUT bspt.f DEPENDS bspt.F  COMMAND cpp -traditional -P bspt.F > bspt.f)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES bspt.f)

# what are the dependencies for display_version.f90, or shuold it build every time?
add_custom_command (OUTPUT display_version.f90 COMMAND ../../SCRIPTS/all/display_version.sh > display_version.f90)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES display_version.f90)

# what are the dependencies for main.f, or shuold it build every time?
add_custom_command (OUTPUT main.f COMMAND cpp -traditional -P main.F > main.f)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES main.f)

# what are the dependencies for porfuncs.f90, or shuold it build every time?
add_custom_command (OUTPUT porfuncs.f90 COMMAND ./porfuncs.csh > porfuncs.f90)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES porfuncs.f90)

################ end More elegany solution needed! ####################

# improve dependencies here and make it build automatically
add_library(mylapack STATIC IMPORTED)
set_property(TARGET mylapack PROPERTY IMPORTED_LOCATION ../../LAPACK/libmylapack.a)

# improve dependencies here and make it build automatically
add_library(myblas STATIC IMPORTED)
set_property(TARGET myblas PROPERTY IMPORTED_LOCATION ../../BLAS/libmyblas.a)

add_subdirectory(AMH)

### now the executables

add_executable(GMIN ${GMIN_SOURCES} main.f mc.f bspt.f display_version.f90 porfuncs.f90 ${DUMMY_CHARMM} ${DUMMY_AMBER} ${DUMMY_AMBER9})
target_link_libraries(GMIN mylapack myblas gfortran)

add_executable(AMHGMIN ${GMIN_SOURCES} main.f mc.f bspt.f display_version.f90 porfuncs.f90 ${DUMMY_CHARMM} ${DUMMY_AMBER9})
target_link_libraries(AMHGMIN amber mylapack myblas gfortran)


  