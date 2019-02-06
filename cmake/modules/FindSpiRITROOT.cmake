# - Find SpiRITROOT instalation
# This module tries to find the SpiRITROOT installation on your system.
#
# Variables defined by this module:
#
#   SpiRITROOT_FOUND               System has SpiRITROOT
#   SpiRITROOT_INCLUDE_DIR         SpiRITROOT include directories: not cached
#   SpiRITROOT_LIBRARY_DIR         The path to where the SpiRITROOT library files are.
#

Message(STATUS "Looking for SpiRITROOT...")


Set(SpiRITROOT_LIBRARY_SEARCHPATH
  ${SPIRITROOTPATH}/build/lib
)

Message(STATUS "Finding -> ${SpiRIT_MODULES}")

Set(SpiRITROOT_FOUND FALSE)

Find_Library(SpiRITROOT_LIBRARY NAMES ${SpiRIT_MODULES}
             PATHS ${SpiRITROOT_LIBRARY_SEARCHPATH}
             NO_DEFAULT_PATH
            )

If(SpiRITROOT_LIBRARY)

  MESSAGE(STATUS "Looking for SpiRITROOT... - found ${SPIRITROOTPATH}/build/lib")

  Set(SpiRITROOT_LIBRARY_DIR ${SPIRITROOTPATH}/build/lib)
  Set(SpiRITROOT_LDFLAGS "-L${SPIRITROOTPATH}/bulid/lib -lSTFormat -lSTGlobal -lSTAna -lSTUtil")

  MESSAGE(STATUS "LDFLAGS ${SpiRITROOT_LDFLAGS}")

  Set(SpiRITROOT_INCLUDE_DIR ${SPIRITROOTPATH}/include)

  MESSAGE(STATUS "includes are in ... ${SpiRITROOT_INCLUDE_DIR}")

  Mark_As_Advanced(SpiRITROOT_LIBRARY_DIR ${SpiRITROOT_INCLUDE_DIR})

  Set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${SpiRITROOT_LIBRARY_DIR})

  Set(SpiRITROOT_FOUND TRUE)

Else(SpiRITROOT_LIBRARY)

  If(SpiRITROOT_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Looking for SpiRITROOT... - Not found!")
  EndIf(SpiRITROOT_FIND_REQUIRED)

EndIf(SpiRITROOT_LIBRARY)
