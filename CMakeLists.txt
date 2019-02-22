cmake_minimum_required(VERSION 3.12)

# Turn on extra compiler warnings:
add_compile_options(-Wall)

# Create a list of the required linked libraries:
set(LINK_LIBS "")

# We need the pow() function, so check that it exists, and link to libm if not:
include(CheckSymbolExists)
check_symbol_exists(pow "" HAS_POW)
if(NOT HAS_POW)
    message("pow does not exist, so linking to libm")
    list(APPEND LINK_LIBS m)
endif()

# Find the HTSlib package:
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(HTSlib REQUIRED)

# Build the bamwindow executable:
add_executable(bamwindow bamwindow.c)
list(APPEND LINK_LIBS ${HTSlib_LIBRARIES})
include_directories(${HTSlib_INCLUDE_DIRS})
target_link_libraries(bamwindow ${LINK_LIBS})

# Install the bamwindow executable:
install(TARGETS bamwindow DESTINATION bin)
