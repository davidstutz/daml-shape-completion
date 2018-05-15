#================================================================
# Install GFlags.
#================================================================

wget -O gflags-v2.2.1.zip https://github.com/gflags/gflags/archive/v2.2.1.zip
unzip gflags-v2.2.1.zip
cd gflags-2.2.1
mkdir build
cd build
cmake ..
make
make DESTDIR=$WORK/dev-box install
cd ..
cd ..