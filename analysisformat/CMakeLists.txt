# Create a library called "libSTAnalysisFormat" which includes the source files given in the array .

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR}
${EIGEN3_INCLUDE_DIR}
${GENFIT2_INCLUDE_DIR}
${BASE_INCLUDE_DIRECTORIES}
# Add all the folders having header files you used in the code.
${CMAKE_SOURCE_DIR}/tasks
${CMAKE_SOURCE_DIR}/analysisformat
${CMAKE_SOURCE_DIR}/flowformat
${SPIRITROOTPATH}/format
${SPIRITROOTPATH}/reco
${SPIRITROOTPATH}/reco/Task
${SPIRITROOTPATH}/reco/Genfit
${SPIRITROOTPATH}/reco/PID
${SPIRITROOTPATH}/reco/Genfit
${SPIRITROOTPATH}/reco/BDCMagneticProjection
${SPIRITROOTPATH}/util
${SPIRITROOTPATH}/global
${SPIRITROOTPATH}/mc
${SPIRITROOTPATH}/field
${SPIRITROOTPATH}/geant4
${SPIRITROOTPATH}/par
${SPIRITROOTPATH}/ana
${MASSFITTER}/
${FAIRROOTPATH}/include
)

include_directories(${INCLUDE_DIRECTORIES})
MESSAGE(STATUS "-->> ${INCLUDE_DIRECTORIES}")

include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${GENFIT2_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
${SPIRITROOTPATH}/build/lib
${CMAKE_SOURCE_DIR}/analysisformat
${CMAKE_SOURCE_DIR}/flowformat
) 

link_directories(${LINK_DIRECTORIES})

set(SRCS
# Add all the source files below this line. Those must have cc for their extension.
STKatana.cc
STParticle.cc
STKParticle.cc
STTriggerArray.cc
STNeuLANDCluster.cc
STNeuLANDHit.cc
STBetheBlochFittingFunction.cc
STBDC.cc
KatanaRoot.cc
STBBFunction.cc
STMassFunction.cc
#STMassCalculator.cc
STMassCalSimpleBB.cc
)

CHANGE_FILE_EXTENSION(*.cc *.hh HEADERS "${SRCS}")

set(HEADERS ${HEADERS})
Set(LINKDEF STAnalysisFormatLinkDef.h)
Set(LIBRARY_NAME STAnalysisFormat)

Set(DEPENDENCIES
  Base
  Stack
# Add all library names your code depending on.
# Basically those are the libraries in the header files' directory
  STFormat
)

GENERATE_LIBRARY()
