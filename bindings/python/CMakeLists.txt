cmake_minimum_required(VERSION 3.8)

project(CONSTITUTIVE_PYTHON VERSION 1.0.0)

# Include python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

cmake_policy(SET CMP0078 OLD)
cmake_policy(SET CMP0086 NEW)

# Include SWIG
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

if(PYTHONLIBS_VERSION_STRING MATCHES "^2.*$" )
    set(CMAKE_SWIG_FLAGS -classic)
else()
    set(CMAKE_SWIG_FLAGS -py3)
endif()

set_source_files_properties(../constitutive.i PROPERTIES CPLUSPLUS ON)

include_directories("${UTILITY_INC_DIR}")
include_directories("${SIMPLEMATH_INC_DIR}")
include_directories(${CONSTITUTIVE_INCLUDE})
include_directories("${PYTHON_INCLUDE_PATH}/../Lib/site-packages/numpy/core/include")

# Add swig module
swig_add_library(constitutive
	LANGUAGE python
	SOURCES ../constitutive.i)
swig_link_libraries(constitutive SimpleMath::SimpleMath)
swig_link_libraries(constitutive Constitutive::Constitutive)
swig_link_libraries(constitutive ${PYTHON_LIBRARIES})

# Bring any examples over to build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_Turon06CZ.py
               ${CMAKE_CURRENT_BINARY_DIR}/test_Turon06CZ.py COPYONLY)


# Files to install with Python
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/constitutive.py
        ${CMAKE_CURRENT_BINARY_DIR}/test_Turon06CZ.py
        DESTINATION python)

# Install the python bindings target
install(
    TARGETS _constitutive
    LIBRARY DESTINATION python
    ARCHIVE DESTINATION python
)

# Declare install target for python
#install(TARGETS swig_example
#        COMMAND "${PYTHON_EXECUTABLE} setup.py"
#        COMPONENT swig-python)