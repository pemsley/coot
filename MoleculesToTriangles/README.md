# MoleculesToTriangles
Routines to triangulate Molecular representations

Linux
cd ClipperAndDependencies/mmdb--2.0.1
mkdir build
cd build
../configure --enable-shared
make -j 8
sudo make install

cd ClipperAndDependencies/fftw-2.1.5
mkdir build
cd build
../configure --enable-shared --with-openmp --enable-float
make -j 8
sudo make install

../configure --enable-ccp4 --enable-shared

cd ClipperAndDependencies/clipper
mkdir build
cd build
CPPFLAGS=-I/usr/local/ccp4/ccp4-7.0/include ../configure --enable-ccp4 --e
nable-shared
make -j 8
sudo make install
