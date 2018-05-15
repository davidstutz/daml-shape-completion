#================================================================
# Install VTK7.
#================================================================

wget http://www.vtk.org/files/release/7.1/VTK-7.1.1.zip
unzip VTK-7.1.1.zip
cd VTK-7.1.1
mkdir build
cd build
cmake ..
make -j4
make DESTDIR=$WORK/dev-box install