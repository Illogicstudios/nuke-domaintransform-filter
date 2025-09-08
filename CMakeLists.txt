cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(NDKBlur)

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

find_package(Nuke REQUIRED)

function(add_custom_plugin PLUGIN_NAME)
    add_nuke_plugin(${PLUGIN_NAME} ${ARGN})
    target_compile_definitions(${PLUGIN_NAME} PRIVATE FN_EXAMPLE_PLUGIN)
endfunction()

add_custom_plugin(DTFilter DTFilter.cpp)
