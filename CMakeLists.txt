cmake_minimum_required(VERSION 3.5)

project(Thunder VERSION 0.1 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_BUILD_TYPE "Debug")

find_package(QT NAMES Qt6 Qt5 COMPONENTS Xml Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Xml Widgets REQUIRED)

# VTK
find_package(VTK REQUIRED)
vtk_module_config(VTK
  vtkFiltersGeometry
  vtkFiltersModeling
  vtkCommonCore
  vtkFiltersSources
  vtkInteractionStyle
  vtkRenderingOpenGL2
  vtkIOGeometry
  vtkGUISupportQt
  vtkRenderingCore
  vtkInteractionWidgets
)
include(${VTK_USE_FILE})

include_directories(
        ./
        # /home/mzy/vtk/VTK-8.2.0/Common/Core/vtkEventData.h
)

set(PROJECT_SOURCES
        main.cpp
        # vtkAxisWidget.cpp
        vtkCustomWidget.cpp
        vtkCustomRepresentation.cpp
        vtkBoxWidget2.cxx
        vtkBoxRepresentation.cxx
        )
include_directories("~/Qt5.12.12/5.12.12/gcc_64/")

add_executable(Thunder
        ${PROJECT_SOURCES})

target_link_libraries(Thunder PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Xml ${VTK_LIBRARIES})

set_target_properties(Thunder PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER my.example.com
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)
