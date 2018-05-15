#================================================================
# Install GLog.
#================================================================

wget -O glog-v0.3.5.zip https://github.com/google/glog/archive/v0.3.5.zip
unzip glog-v0.3.5.zip
cd glog-0.3.5
mkdir build
cd build
cmake ..
make
make DESTDIR=$WORK/dev-box install
cd ..
cd ..