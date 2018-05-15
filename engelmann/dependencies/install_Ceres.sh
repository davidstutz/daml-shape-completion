#================================================================
# Install Ceres.
#================================================================

echo "Dependency: SuiteSparse (and its dependecies, use install_SuiteSparse.sh)."
echo "Dependency: GLog (use install_GLog.sh)."
echo "Dependency: GFlags (use install_GFlags.sh)."
echo "Dependency: ATLAS (shoul dbe pre-installed globally, check /usr/lib/libatlas.so or locate libatlas)."
echo "Dependency: Eigen3 (use install_Eigen.sh)."
echo "Also make sure that the hint paths for CMake are set correctly (in this file);"
echo "The above dependecies are assumed to be installed in $WORK/dev-box/usr/local/!"
read -p "Press ENTER after making sure that all dependencies are met."

wget http://ceres-solver.org/ceres-solver-1.13.0.tar.gz 
tar -xvzf ceres-solver-1.13.0.tar.gz 
cd ceres-solver-1.13.0
mkdir build
cmake \
-DGFLAGS_INCLUDE_DIR_HINTS="$WORK/dev-box/usr/local/include/" \
-DGFLAGS_LIBRARY_DIR_HINTS="$WORK/dev-box/usr/local/lib/" \
-DGLOG_INCLUDE_DIR_HINTS="$WORK/dev-box/usr/local/include/" \
-DGLOG_LIBRARY_DIR_HINTS="$WORK/dev-box/usr/local/lib/" \
-DSUITESPARSE_INCLUDE_DIR_HINTS="$WORK/dev-box/usr/local/include/" \
-DSUITESPARSE_LIBRARY_DIR_HINTS="$WORK/dev-box/usr/local/lib/" \
-DEIGEN_INCLUDE_DIR_HINTS="$WORK/dev-box/usr/local/include/" \
..
make -j4
make DESTDIR=$WORK/dev-box install
cd ..
cd ..