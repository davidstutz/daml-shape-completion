#================================================================
# Install OpenCV.
#================================================================

echo "Remember to make sure to not have more than one distinct version of OpenCV installed."
echo "If you installed a different version previously, use make uninstall DESTDIR=$WORK/dev-box install to uninstall (within the corresponding build directory)."
read -p "Press ENTER to continue."

wget -O opencv-2.4.13.3.zip https://github.com/opencv/opencv/archive/2.4.13.3.zip
unzip opencv-2.4.13.3.zip
cd opencv-2.4.13.3
mkdir build
cd build
cmake ..
make -j4
make DESTDIR=$WORK/dev-box install
cd ..
cd ..