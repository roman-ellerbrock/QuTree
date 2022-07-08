
# Set up BLAS / LAPACK
# BLAS++ & LAPACK++ per default install to /opt/slate/lib, thus added a hint
find_library(BLAS_LIBRARIES
	NAMES blaspp
	HINTS /opt/slate/lib 
	REQUIRED)
find_library(LAPACK_LIBRARIES
	NAMES lapackpp
	HINTS /opt/slate/lib 
	REQUIRED)
message("BLAS++ libraries: ${BLAS_LIBRARIES}")
message("LAPACK++ libraries: ${LAPACK_LIBRARIES}")

