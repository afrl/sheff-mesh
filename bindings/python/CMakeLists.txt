cmake_minimum_required(VERSION 3.8)

project(MESH_PYTHON VERSION 1.0.0)

# Include python

# This is a hack to fix how CMake detects Python.  CMake just looks at the location of the
# "python" command, but this may be a shim (see pyenv).  Instead, the location is queried
# from python itself before CMake tries to find Python, forcing it to use this location,
# instead of the shim.
execute_process(COMMAND cmd;/c;python -c "import sys;print(sys.executable)"
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE python_loc_result
                OUTPUT_VARIABLE PYTHON_EXECUTABLE)
message(STATUS "Python command result: ${python_loc_result}")
message(STATUS "Python executable location: ${PYTHON_EXECUTABLE}")

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

# Include SWIG
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

if(PYTHONLIBS_VERSION_STRING MATCHES "^2.*$" )
    set(CMAKE_SWIG_FLAGS -classic)
else()
    set(CMAKE_SWIG_FLAGS -py3)
endif()

set_source_files_properties(../mesh.i PROPERTIES CPLUSPLUS ON)

include_directories("${UTILITY_INC_DIR}")
include_directories("${SIMPLEMATH_INC_DIR}")
include_directories("${MESH_INCLUDE}")

execute_process(COMMAND cmd;/c;python -c "import numpy;print(numpy.__file__)"
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE numpy_loc_result
                OUTPUT_VARIABLE NUMPY_LOCATION)
message(STATUS "Numpy command result: ${numpy_loc_result}")
message(STATUS "Numpy location: ${NUMPY_LOCATION}")
get_filename_component(NUMPY_DIR "${NUMPY_LOCATION}" DIRECTORY)
message(STATUS "Numpy directory: ${NUMPY_DIR}")

include_directories("${NUMPY_DIR}/core/include")

# Add swig module
swig_add_library(mesh
	LANGUAGE python
	SOURCES ../mesh.i)
swig_link_libraries(mesh SimpleMath::SimpleMath)
swig_link_libraries(mesh Mesh::Mesh)
swig_link_libraries(mesh ${PYTHON_LIBRARIES})

# Bring any examples over to build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_Mesh.py
               ${CMAKE_CURRENT_BINARY_DIR}/test_Mesh.py COPYONLY)


# Files to install with Python
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/mesh.py
        ${CMAKE_CURRENT_BINARY_DIR}/test_Mesh.py
        DESTINATION python)

# Install the python bindings target
install(
    TARGETS _mesh
    LIBRARY DESTINATION python
    ARCHIVE DESTINATION python
)

# Declare install target for python
#install(TARGETS swig_example
#        COMMAND "${PYTHON_EXECUTABLE} setup.py"
#        COMPONENT swig-python)