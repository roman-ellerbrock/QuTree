# How to install QuTree on Fire CPU nodes

#-----Load required modules----#
# (if not already done)
ml CMake/3.7.1-intel-2016.0.109
ml GCCcore/9.3.0
ml intel/2017.8.262


#------Set up install dir------#
# (if it does not already exist)
mkdir [PATH_TO_INSTALL_DIR]         # e.g. /data/kgjohn/usr
mkdir [PATH_TO_INSTALL_DIR]/include # e.g. /data/kgjohn/usr/include


#---------Install QuTree-------#
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=[PATH_TO_INSTALL_DIR]
make install
