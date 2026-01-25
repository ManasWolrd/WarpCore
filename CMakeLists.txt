cmake_minimum_required(VERSION 3.22.1)

# ----------------------------------------
# 插件信息，改了记得改src/global.hpp
# ----------------------------------------
set(TEAM_NAME LModel)
set(PLUGIN_REPO_URL "https://github.com/team_name/repo_name")
set(PLUGIN_NAME template)
set(PLUGIN_VERSION 0.0.1)
set(PLUGIN_BUNDLE_ID "com.l_model.${PLUGIN_NAME}")

option(BUILD_STANDALONE "Build Standalone plugin format" ON)
option(BUILD_VST3 "Build VST3 plugin format" ON)
option(BUILD_LV2 "Build LV2 plugin format" ON)

# cmake设置

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>" CACHE STRING "Runtime Library")
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "Build for 10.13")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(${PLUGIN_NAME} VERSION ${PLUGIN_VERSION})

add_compile_definitions(PROJECT_VERSION="${PROJECT_VERSION}")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_XCODE_GENERATE_SCHEME OFF)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

string(SUBSTRING plugin_codes ${PLUGIN_NAME} 0 4)
string(SUBSTRING plugin_manufacturer_code ${TEAM_NAME} 0 4)

add_subdirectory(libs/JUCE)

set(plugin_formats AU)

if(BUILD_STANDALONE)
    list(APPEND plugin_formats Standalone)
endif()

if(BUILD_VST3)
    list(APPEND plugin_formats VST3)
endif()

if (BUILD_LV2 AND UNIX AND NOT APPLE)
    list(APPEND plugin_formats LV2)
endif()

# 插件信息在这里也可以设置
juce_add_plugin(${PROJECT_NAME}
    COMPANY_NAME ${TEAM_NAME}
    IS_SYNTH FALSE
    NEEDS_MIDI_INPUT False
    NEEDS_MIDI_OUTPUT false
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    JUCE_VST3_CAN_REPLACE_VST2 FALSE
    COPY_PLUGIN_AFTER_BUILD TRUE
    PLUGIN_MANUFACTURER_CODE ${plugin_manufacturer_code}
    PLUGIN_CODE ${plugin_codes}
    FORMATS ${plugin_formats}
    PRODUCT_NAME ${PLUGIN_NAME}
    BUNDLE_ID ${PLUGIN_BUNDLE_ID}
    LV2URI ${PLUGIN_REPO_URL}
)

# 添加源文件
file(GLOB_RECURSE src
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
target_sources(${PROJECT_NAME} PRIVATE ${src})

# 头文件目录
target_include_directories(${PROJECT_NAME}
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

# juce编译选项
target_compile_definitions(${PROJECT_NAME}
    PUBLIC
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
)

# 要用 #include <juce.h> 去掉这个注释
# juce_generate_juce_header(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        juce::juce_dsp
        juce::juce_core
        juce::juce_graphics
        juce::juce_gui_basics
        juce::juce_audio_utils
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags
)

# 添加资源
file(GLOB_RECURSE res "${CMAKE_CURRENT_SOURCE_DIR}/resources/*.*")
if (NOT "${res}" STREQUAL "")
    juce_add_binary_data(${PROJECT_NAME}_res SOURCES ${res})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_res)
endif()

if(APPLE)
    target_compile_definitions(${PROJECT_NAME} PUBLIC JUCE_AU=1)
endif()
