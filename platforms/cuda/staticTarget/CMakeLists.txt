#
# Include CUDA related files.
#
INCLUDE(${FINDCUDA_DIR}/FindCuda.cmake)
INCLUDE_DIRECTORIES(${CUDA_INCLUDE})
LINK_DIRECTORIES(${CUDA_TARGET_LINK})
FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB src_files ${CMAKE_SOURCE_DIR}/platforms/cuda/${subdir}/src/*.cu  ${CMAKE_SOURCE_DIR}/platforms/cuda/src/*/*.cu)
    FOREACH(file ${src_files})
        FILE(RELATIVE_PATH file ${CMAKE_SOURCE_DIR}/platforms/cuda ${file})
        SET(SOURCE_FILES         ${SOURCE_FILES}         ${file})   #append
    ENDFOREACH(file)
    CUDA_INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/cuda/../${subdir}/include)
ENDFOREACH(subdir)

CUDA_INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/cuda/../src)

CUDA_ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${STATIC_TARGET} debug ${OPENMM_LIBRARY_NAME}_static_d optimized ${OPENMM_LIBRARY_NAME}_static)
