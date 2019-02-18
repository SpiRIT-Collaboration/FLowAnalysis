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
${SPIRITROOTPATH}/include
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
${CMAKE_SOURCE_DIR}/tasks
) 

link_directories(${LINK_DIRECTORIES})
MESSAGE(STATUS "SpiRITAnalysis LINK_DIRECTORIES -> ${LINK_DIRECTORIES} <--")

set(SRCS
# Add all the source files below this line. Those must have cc for their extension.
STBigRIPSTask.cc
STSpiRITTPCTask.cc
)

CHANGE_FILE_EXTENSION(*.cc *.hh HEADERS "${SRCS}")

set(HEADERS ${HEADERS})
Set(LINKDEF STAnalysisTaskLinkDef.h)
Set(LIBRARY_NAME STAnalysisTasks)
Set(DEPENDENCIES
  Base
  Stack
# Add all library names your code depending on.
# Basically those are the libraries in the header files' directory
  STGlobal
  STAna
  STAnalysisFormat
  STFlowAnalysis
)

#target_link_libraries(${LIBRARY_NAME} ${SpiRIT_MODULES})


GENERATE_LIBRARY()
