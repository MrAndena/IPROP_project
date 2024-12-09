# ##
# #  CMake script for the time-dependent Navier-Stokes solver:
# ##

# # Set the name of the project and target:
# SET(TARGET "fluid")

# # Declare all source files the target consists of. Here, this is only
# # the one step-X.cc file, but as you expand your project you may wish
# # to add other source files as well. If your project becomes much larger,
# # you may want to either replace the following statement by something like
# #  FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
# #  FILE(GLOB_RECURSE TARGET_INC  "include/*.h")
# #  SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# # or switch altogether to the large project CMakeLists.txt file discussed
# # in the "CMake in user projects" page accessible from the "User info"
# # page of the documentation.
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -g")
# #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -v -Wl, --warn-unresolved-symbols")
# FILE(GLOB_RECURSE TARGET_SRC "time_dependent_navier_stokes.cc")
# FILE(GLOB_RECURSE TARGET_INC  "MyDataStruct.hpp" "Time.hpp" "BoundaryValues.hpp" "json.hpp" "BlockSchurPreconditioner.hpp" "BlockSchurPreconditioner_impl.hpp" "InsIMEX.hpp" "InsIMEX_impl.hpp" )
# SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})

# #SET(TARGET_SRC
#  # MyDataStruct.hpp 
#   #Time.hpp
#   #BoundaryValues.hpp
#  # json.hpp

#   #BlockSchurPreconditioner.hpp 
#   #BlockSchurPreconditioner.cpp

#   #InsIMEX.hpp 
#   #InsIMEX.cpp

#   #time_dependent_navier_stokes.cc
#   #)

# # Usually, you will not need to modify anything beyond this point...

# CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)   

# #FIND_PACKAGE(deal.II 9.5
# #  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
# #  )

# FIND_PACKAGE(deal.II 9.5.1
#     # HINTS /uu/sw/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.3.1
#     HINTS /opt/mox/mk/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.5.1
#     #HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
#     )

# IF(NOT ${deal.II_FOUND})
#   MESSAGE(FATAL_ERROR "\n"
#     "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
#     "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
#     "or set an environment variable \"DEAL_II_DIR\" that contains this path."
#     )
# ENDIF()

# # FIND_PACKAGE(MPI REQUIRED)
# # SET (CMAKE_C_COMPILER ${MPI_C_COMPILER})
# # SET (CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# # # Specifica il numero di processori
# # SET(MPIEXEC_NUMPROC_FLAG "-n")
# # SET(MPIEXEC_MAX_NUMPROCS "1")  # Cambia il numero di processori a tuo piacimento

# # # Aggiungi gli flags MPI
# # SET(MPIEXEC_PREFLAGS "${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}")

# #
# # Are all dependencies fulfilled?
# #
# IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_P4EST OR DEAL_II_PETSC_WITH_COMPLEX) # keep in one line
#   MESSAGE(FATAL_ERROR "
# Error! This tutorial requires a deal.II library that was configured with the following options:
#     DEAL_II_WITH_PETSC = ON
#     DEAL_II_PETSC_WITH_COMPLEX = OFF
#     DEAL_II_WITH_P4EST = ON
# However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
#     DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
#     DEAL_II_PETSC_WITH_COMPLEX = ${DEAL_II_PETSC_WITH_COMPLEX}
#     DEAL_II_WITH_P4EST = ${DEAL_II_WITH_P4EST}
# which conflict with the requirements."
#     )
# ENDIF()


# DEAL_II_INITIALIZE_CACHED_VARIABLES()
# PROJECT(${TARGET})
# DEAL_II_INVOKE_AUTOPILOT()


# #Per eseguire mettersi nella cartelle build e fare:
# # a) cmake .. (solo se si ha modificato il file CMakeLists.txt, se no partire dal punto b)
# # b) cmake --build .  oppure   make 
# # c) ./fluid   oppure    mpiexec -n 4 ./fluid



# CMake minimum version
cmake_minimum_required(VERSION 3.13.4)

# Set the name of the project and target:
set(TARGET "fluid")

# Deal.ii
FIND_PACKAGE(deal.II 9.5.1
    HINTS /opt/mox/mk/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.5.1
)

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
  )
ENDIF()

# Ensure deal.II is correctly configured
IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_P4EST OR DEAL_II_PETSC_WITH_COMPLEX)
  MESSAGE(FATAL_ERROR "
Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_PETSC = ON
    DEAL_II_PETSC_WITH_COMPLEX = OFF
    DEAL_II_WITH_P4EST = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options:
    DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
    DEAL_II_PETSC_WITH_COMPLEX = ${DEAL_II_PETSC_WITH_COMPLEX}
    DEAL_II_WITH_P4EST = ${DEAL_II_WITH_P4EST}
which conflict with the requirements."
  )
ENDIF()

# Add include directories for custom headers
include_directories(${CMAKE_SOURCE_DIR}/include)

# Declare source files
FILE(GLOB_RECURSE TARGET_SRC "src/*.cc" "include/*.hpp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -g")

# Initialize deal.II variables and project
DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})

# This macro auto-configures the build for deal.II projects
DEAL_II_INVOKE_AUTOPILOT()

# NEW command
# make custom_clean per eliminare i file specificati
add_custom_target(custom_clean
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/compile_commands.json     #delete commands.json
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/output/            #delete all simulation results
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/output/              #create new output directory
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}                    #delete build direcotry
    COMMENT "Custom clean: Deleting compile_commands.json, build directory and the numerical simulations"
)

# make clean_output per eliminare solo gli output
add_custom_target(clean_output
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/output/            #delete all simulation results
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/output/              #create new output directory
    COMMENT "clean output: Deleting all the numerical simulations"
)