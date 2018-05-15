 
wget -O openblas-v0.2.20.zip http://github.com/xianyi/OpenBLAS/archive/v0.2.20.zip
unzip openblas-v0.2.20.zip
cd OpenBLAS-0.2.20
mkdir build
cd build
cmake ..
make -j4
make DESTDIR=$WORK/dev-box install
cd ..
cd ..