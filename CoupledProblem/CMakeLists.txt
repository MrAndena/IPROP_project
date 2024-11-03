# CMake minimum version
cmake_minimum_required(VERSION 3.13.4)

# Set the name of the project and target:
set(TARGET "MyExecutable")

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
FILE(GLOB_RECURSE TARGET_SRC "src/*.cpp" "include/*.hpp")
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





