# The name of our project is “MYPROJ”.  CMakeLists files in this project can
# refer to the root source directory of the project as ${MYPROJ_SOURCE_DIR} 
# or as ${CMAKE_SOURCE_DIR} and to the root binary directory of the project as 
# ${MYPROJ_BINARY_DIR} or ${CMAKE_BINARY_DIR}.
# This difference is important for the base classes which are in FAIRROOT 
# and the experiment part.

# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2 FATAL_ERROR)

### CMP0025   Compiler id for Apple Clang is now AppleClang.
### CMP0042   MACOSX_RPATH is enabled by default.

foreach(p
  CMP0025 # CMake 3.0
  CMP0042 # CMake 3.0
  CMP0054 # CMake 3.9
  )
  if(POLICY ${p})
  cmake_policy(SET ${p} NEW)
  endif()
endforeach()

# Set name of our project to “MYPROJ". Has to be done
# after check of cmake version since this is a new feature                  
project(SpiRITAnalysis)

#In case you need Fortran 
#ENABLE_LANGUAGE(Fortran)

# Check for needed environment variables
IF(NOT DEFINED ENV{FAIRROOTPATH})
  MESSAGE(FATAL_ERROR "You did not define the environment variable FAIRROOTPATH which is needed to find FairRoot. Please set this variable and execute cmake again.")
ENDIF(NOT DEFINED ENV{FAIRROOTPATH})

IF(NOT DEFINED ENV{SIMPATH})
   MESSAGE(FATAL_ERROR "You did not define the environment variable SIMPATH which is nedded to find the external packages. Please set this variable and execute cmake again.") 
ENDIF(NOT DEFINED ENV{SIMPATH})

IF(NOT DEFINED ENV{SPIRITROOTPATH})
   MESSAGE(FATAL_ERROR "You did not define the environment variable SPIRITROOTPATH which is nedded to find the external packages. Please set this variable and execute cmake again.") 
ENDIF(NOT DEFINED ENV{SPIRITROOTPATH})

SET(SIMPATH        $ENV{SIMPATH})
SET(FAIRROOTPATH   $ENV{FAIRROOTPATH})
SET(SPIRITROOTPATH $ENV{SPIRITROOTPATH})
SET(SPIRITANAPATH  $ENV{SPIRITANAPATH})

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ 
# is checked
set(CMAKE_MODULE_PATH "${SPIRITROOTPATH}/cmake/modules")
set(CMAKE_MODULE_PATH "${FAIRROOTPATH}/share/fairbase/cmake/modules"  ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${FAIRROOTPATH}/share/fairbase/cmake/modules_old"  ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

execute_process(COMMAND echo ${CMAKE_MODULE_PATH})

Set(CheckSrcDir "${FAIRROOTPATH}/share/fairbase/cmake/checks")
find_package(FairRoot)



# Load some basic macros which are needed later on
include(FairMacros)
include(WriteConfigFile)
include(CTest)
include(CheckCompiler)
#include(CheckFortran)

#Check the compiler and set the compile and link flags
If(NOT CMAKE_BUILD_TYPE)
  Message("Set BuildType DEBUG")
  set(CMAKE_BUILD_TYPE Debug)
EndIf(NOT CMAKE_BUILD_TYPE)
Check_Compiler()


set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(INCLUDE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/include")
Set(VMCWORKDIR ${CMAKE_SOURCE_DIR})
Option(USE_PATH_INFO "Information from PATH and LD_LIBRARY_PATH are used." ON)
If(USE_PATH_INFO)
  Set(PATH $ENV{PATH})
  If (APPLE)
    Set(LD_LIBRARY_PATH $ENV{DYLD_LIBRARY_PATH})
  Else (APPLE)
    Set(LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})
    MESSAGE(" NOT apple " ${LD_LIBRARY_PATH} )
  EndIf (APPLE)
EndIf(USE_PATH_INFO)


# Check if the user wants to build the project in the source
# directory
CHECK_OUT_OF_SOURCE_BUILD()

# Check if we are on an UNIX system. If not stop with an error
# message
IF(NOT UNIX)
  MESSAGE(FATAL_ERROR "You're not on an UNIX system. The project was up to now only tested on UNIX systems, so we break here. If you want to go on please edit the CMakeLists.txt in the source directory.")
ENDIF(NOT UNIX)  

# Check if the external packages are installed into a separate install
# directory
CHECK_EXTERNAL_PACKAGE_INSTALL_DIR()

# searches for needed packages
# REQUIRED means that cmake will stop if this packages are not found
# For example the framework can run without GEANT4, but ROOT is
# mandatory

find_package(ROOT 5.32.00 REQUIRED)
find_package(GEANT3 REQUIRED)
find_package(GEANT4 REQUIRED)
find_package(GEANT4VMC)
find_package(CLHEP)
find_package(Eigen3 REQUIRED)
find_package(GENFIT2 REQUIRED)
find_package(ANAROOT REQUIRED)
find_package(RAVE REQUIRED)
find_package(OpenMP)
find_package(PythonInterp 2.7 REQUIRED)

If(NOT APPLE)
  set(PYTHON_LIBRARY ${PYTHON_EXECUTABLE})
  set(PYTHON_INCLUDE_DIR ${PYTHON_EXECUTABLE})
  string(REGEX REPLACE "bin(.)+" "lib64" PYTHON_LIBRARY ${PYTHON_LIBRARY})
  string(REGEX REPLACE "bin(.)+" "include/python${PYTHON_VERSION_STRING}" PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR})
EndIf(NOT APPLE)
find_package(PythonLibs 2.7 REQUIRED)

Set(Boost_NO_SYSTEM_PATHS TRUE)
Set(Boost_NO_BOOST_CMAKE TRUE)
If(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
  set(BOOST_ROOT ${SIMPATH})
  set(GSL_DIR ${SIMPATH})
Else(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
  set(BOOST_ROOT ${SIMPATH}/basics/boost)
  set(GSL_DIR ${SIMPATH}/basics/gsl)
EndIf(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
Message("-- Looking for Boost ...")
# If an older version of boost is found both of the variables below are
# cached and in a second cmake run, a good boost version is found even 
# if the version is to old. 
# To overcome this problem both variables are cleared before checking
# for boost.
Unset(Boost_INCLUDE_DIR CACHE)
Unset(Boost_LIBRARY_DIRS CACHE)
find_package(Boost 1.41 COMPONENTS thread system timer program_options random)
If (Boost_FOUND)
  Set(Boost_Avail 1)
  MESSAGE( "bf  -> " ${LD_LIBRARY_PATH} )
  MESSAGE( "bst -> " ${Boost_LIBRARY_DIR} )
  Set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${Boost_LIBRARY_DIR})
  MESSAGE( "<--:" ${LD_LIBRARY_PATH} )
Else (Boost_FOUND)
  Set(Boost_Avail 0)
EndIf (Boost_FOUND)

# set a variable which should be used in all CMakeLists.txt
# Defines all basic include directories from fairbase
SetBasicVariables()

# Set the library version in the main CMakeLists.txt
SET(FAIRROOT_MAJOR_VERSION 0)
SET(FAIRROOT_MINOR_VERSION 0)
SET(FAIRROOT_PATCH_VERSION 0)
SET(FAIRROOT_VERSION "${FAIRROOT_MAJOR_VERSION}.${FAIRROOT_MINOR_VERSION}.${FAIRROOT_PATCH_VERSION}")
If(ROOT_FOUND_VERSION LESS 59999)
  SET(FAIRROOT_LIBRARY_PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES}
    VERSION "${FAIRROOT_VERSION}"
    SOVERSION "${FAIRROOT_MAJOR_VERSION}"
    SUFFIX ".so"
  )
Else()
  SET(FAIRROOT_LIBRARY_PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES}
    VERSION "${FAIRROOT_VERSION}"
    SOVERSION "${FAIRROOT_MAJOR_VERSION}"
  )
EndIf()

Set(SpiRIT_MODULES STFormat STGlobal STAna STUtil)
find_package(SpiRITROOT REQUIRED)

Generate_Version_Info()


SET(_LIBDIR ${CMAKE_BINARY_DIR}/lib)
SET(LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})
SET(LD_LIBRARY_PATH  ${_LIBDIR} ${LD_LIBRARY_PATH} ${PYTHON_LIBRARY})

# OpenMP configuration
If(OPENMP_FOUND)
  MESSAGE("-- Compiled with OpenMP")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
EndIf(OPENMP_FOUND)
#install(DIRECTORY geometry DESTINATION pnd_install 
#          )


# Check if the compiler support specific C++11 features
# Up to now this is only a check since the code does not use
# any of the features of the new standard
include(CheckCXX11Features)

IF(HAS_CXX11_SHAREDPOINTER)
  Add_Definitions(-DHAS_SHAREDPOINTER)
ENDIF(HAS_CXX11_SHAREDPOINTER)

# Recurse into the given subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.

add_subdirectory(analysisformat)
add_subdirectory(flowformat)
add_subdirectory(tasks)

WRITE_CONFIG_FILE(config.sh)
 
configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake 
               ${CMAKE_BINARY_DIR}/CTestCustom.cmake
              )

configure_file(${VMCWORKDIR}/cmake/STLog.cmake.in ${VMCWORKDIR}/cmake/STLog.cmake @ONLY)
add_custom_target(STLOG ALL
                  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/STLog.cmake
                  VERBATIM)

file(MAKE_DIRECTORY ${VMCWORKDIR}/macros/data)
