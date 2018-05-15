#================================================================
# Install Eigen3.
#================================================================

# Corresponding cmake file can be found in cmake/.
# Note that in find_path/find_library, NO_CMAKE_SYSTEM_PATH can be used
# to prevent using system paths first in case an old Eigen installation is already
# available.
hg clone https://bitbucket.org/eigen/eigen
cd eigen
mkdir build
cd build
cmake ..
make DESTDIR=$WORK/dev-box install
cd ..
cd ..
