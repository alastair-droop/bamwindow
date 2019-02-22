cmake_minimum_required(VERSION 3.2)

# Turn on extra compiler warnings:
add_compile_options(-Wall)

# Find the HTSlib package:
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(HTSlib REQUIRED)

# Build the bamwindow executable:
add_executable(bamwindow bamwindow.c)
include_directories(${HTSlib_INCLUDE_DIRS})
target_link_libraries(bamwindow ${HTSlib_LIBRARIES})

# Install the bamwindow executable:
install(TARGETS bamwindow DESTINATION bin)
