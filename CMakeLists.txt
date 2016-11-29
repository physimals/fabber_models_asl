add_library(fabber_asl SHARED fwdmodel_asl_rest.cc asl_models.cc)

find_path(FABBER_INCLUDE_DIR fabber_core/fwdmodel.h
            HINTS $ENV{FSLDIR}/include REQUIRED)
find_path(NEWMAT_INCLUDE_DIR newmat.h
            HINTS $ENV{FSLDIR}/extras/include/newmat REQUIRED)
find_path(FSL_INCLUDE_DIR newimage/newimage.h
            HINTS $ENV{FSLDIR}/include REQUIRED)
find_path(FSL_EXT_INCLUDE_DIR libprob/libprob.h
            HINTS $ENV{FSLDIR}/extras/include REQUIRED)
find_path(BOOST_INCLUDE_DIR boost/shared_ptr.hpp
          HINTS $ENV{FSLDIR}/extras/include/boost REQUIRED)

find_library(FABBER_LIBRARY NAMES fabbercore libfabbercore
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(PROB_LIBRARY NAMES prob libprob
             HINTS $ENV{FSLDIR}/extras/lib REQUIRED )
find_library(NEWMAT_LIBRARY NAMES newmat libnewmat
             HINTS $ENV{FSLDIR}/extras/lib REQUIRED )
find_library(MISCMATHS_LIBRARY NAMES miscmaths libmiscmaths
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(UTILS_LIBRARY NAMES utils libutils
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(UTILS_LIBRARY NAMES utils libutils
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(NIFTIIO_LIBRARY NAMES niftiio libniftiio
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(ZNZ_LIBRARY NAMES znz libznz
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(Z_LIBRARY NAMES z libz zlib zlibd 
               HINTS $ENV{FSLDIR}/lib REQUIRED )

include_directories(SYSTEM ${FABBER_INCLUDE_DIR} ${NEWMAT_INCLUDE_DIR} ${FSL_INCLUDE_DIR} ${FSL_EXT_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})
target_link_libraries(fabber_asl ${PROB_LIBRARY} ${NEWMAT_LIBRARY} ${FABBER_LIBRARY} ${MISCMATHS_LIBRARY} ${UTILS_LIBRARY} ${NIFTIIO_LIBRARY} ${ZNZ_LIBRARY} ${Z_LIBRARY})

if (WIN32)
#GENERATE_EXPORT_HEADER( MyLibrary
#             BASE_NAME MyLibrary
#             EXPORT_MACRO_NAME MyLibrary_EXPORT
#             EXPORT_FILE_NAME MyLibrary_Export.h
#             STATIC_DEFINE MyLibrary_BUILT_AS_STATIC)
endif(WIN32)
