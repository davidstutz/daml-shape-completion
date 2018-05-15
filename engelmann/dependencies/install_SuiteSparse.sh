#================================================================
# Install SuiteSparse.
#================================================================

echo "Dependency: OpenBLAS (install_OpenBLAS.sh)."
echo "The above dependecies are assumed to be installed in $WORK/dev-box/usr/local/!"
read -p "Press ENTER after making sure that all dependencies are met."

wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.5.6.tar.gz 
tar -xvzf SuiteSparse-4.5.6.tar.gz 
cd SuiteSparse

CUDA_PATH_='which nvcc 2>/dev/null | sed "s/\/bin\/nvcc//"'
echo "In order to find CUDA, the configuration file of SuiteSparse needs to be adapted."
echo "Go to SuiteSparse/SuiteSparse_config and open SuiteSparse_config.mk."
echo "Search for libcudart.so and replace lib64/libcudart.so by the relevant path."
echo "In particular, the configuration will look for the nvcc binary and get the corresponding path."
echo "Then, CUDA_PATH is set to the corresponding path, so that only the relative path from CUDA_PATH needs to be added."
echo "Currently, CUDA_PATH will be CUDA_PATH=\"$CUDA_PATH_\""
echo "Finally, in the same file, make sure that OpenBLAS is set correctly."
echo "Search for -lopenblas and replace it by the appropriate path of the locally installed OpenBLAS."
echo "This might look as follows: -L/BS/dstutz/work/dev-box/usr/local/lib -lopenblas."
echo "(Or leave it, if OpenBLAS is/can be installed through apt-get.)"
read -p "Press ENTER to continue installation."

make
make install INSTALL=$WORK/dev-box/usr/local/
cd ..