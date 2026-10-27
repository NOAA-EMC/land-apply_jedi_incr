################################################################################
# Dependencies
################################################################################

find_package(OpenMP COMPONENTS C Fortran)
find_package(MPI REQUIRED COMPONENTS C Fortran)
find_package(NetCDF REQUIRED COMPONENTS Fortran )

#add_definitions ( -DCPP=0 )

set( JEDIINCR_LINKER_LANGUAGE Fortran )

# public_libs: snow and soil incr
set( PUBLIC_LIBS
    NetCDF::NetCDF_Fortran
    MPI::MPI_Fortran
    OpenMP::OpenMP_C
    OpenMP::OpenMP_Fortran
)

set( MODULE_DIR ${CMAKE_BINARY_DIR}/module )

################################################################################
# jediincr_soil
################################################################################

list( APPEND sorc_files_soilinc
    ./soil_increments.f90
    ./set_soilveg_snippet.f90
    ./sflx_snippet.f90
)

ecbuild_add_library( TARGET          jediincr_soil
                     SOURCES         ${sorc_files_soilinc}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${JEDIINCR_LINKER_LANGUAGE}
                    )

target_link_libraries(jediincr_soil PUBLIC ${PUBLIC_LIBS})

set_target_properties(jediincr_soil PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR}/soil_inc)
target_include_directories(jediincr_soil PUBLIC $<BUILD_INTERFACE:${MODULE_DIR}/soil_inc>)

set( soil_exe_files ./apply_incr_noahmp_soil.f90 )

ecbuild_add_executable( TARGET          apply_soil_incr.x
                        SOURCES         ${soil_exe_files}
                        LIBS            jediincr_soil
                        LINKER_LANGUAGE ${JEDIINCR_LINKER_LANGUAGE}
                       )
################################################################################
