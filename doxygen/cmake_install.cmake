# Install script for directory: /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdocumentationx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc/doxygen" TYPE FILE MESSAGE_LAZY FILES "/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/deal.tag")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdocumentationx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc/doxygen/deal.II" TYPE FILE MESSAGE_LAZY FILES
    "/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/deal.ico"
    "/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/custom.js"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdocumentationx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc/doxygen" TYPE DIRECTORY MESSAGE_LAZY FILES "/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/deal.II")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial/cmake_install.cmake")
  include("/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/code-gallery/cmake_install.cmake")

endif()

