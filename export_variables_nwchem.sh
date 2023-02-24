#!/usr/bin/sh

#from https://github.com/nwchemgit/nwchem/wiki/Compiling-NWChem
export NWCHEM_TOP=/home/mdi0316/CODES/nwchem-6.8
export NWCHEM_TARGET='LINUX64'
export ARMCI_NETWORK=OPENIB

#$ mpif90 -show
#gfortran -I/cm/shared/apps/openmpi/gcc/64/1.10.7/include -pthread -I/cm/shared/apps/openmpi/gcc/64/1.10.7/lib64 -Wl,-rpath -Wl,/cm/shared/apps/openmpi/gcc/64/1.10.7/lib64 -Wl,--enable-new-dtags -L/cm/shared/apps/openmpi/gcc/64/1.10.7/lib64 -lmpi_usempi -lmpi_mpifh -lmpi

export USE_MPI=y
export LIBMPI="-pthread -rpath -lmpi_usempi -lmpi_mpifh -lmpi"
export MPI_LOC="/cm/shared/apps/openmpi/gcc/64/1.10.7"
export MPI_INCLUDE="$MPI_LOC"/include
export MPI_LIB="$MPI_LOC"/lib
export LD_LIBRARY_PATH=$MPI_LIB:$LD_LIBRARY_PATH


export NWCHEM_MODULES="all"
export USE_NOFSCHECK=TRUE
export USE_NOIO=TRUE
export LIB_DEFINES=-DDFLT_TOT_MEM=16777216
export MRCC_METHODS=TRUE

export PYTHONHOME=/home/mdi0316/anaconda3
export PYTHONVERSION=3.6
export PYTHON_EXE=/home/mdi0316/anaconda3/bin/python
export USE_PYTHONCONFIG=Y
export USE_PYTHON64=Y
export PYTHONLIBTYPE=so
export PYTHONCONFIGDIR=config-x86_64-linux-gnu

#from INSTALL 
export LARGE_FILES=TRUE
export USE_NOFSCHECK=TRUE
#export TCGRSH=/usr/local/bin/ssh

export NWCHEM_EXECUTABLE=$NWCHEM_TOP/bin/LINUX64/nwchem
export NWCHEM_BASIS_LIBRARY=$NWCHEM_TOP/src/basis/libraries/
export NWCHEM_NWPW_LIBRARY=$NWCHEM_TOP/src/nwpw/libraryps/

#export NWCHEM_LONG_PATHS=Y
#export FC=gfortran
#export CC=gcc
#export CXX=g++

export BLASOPT="-L/usr/lib64 -llapack -lblas"

echo "NWCHEM_TOP           = $NWCHEM_TOP"
echo "NWCHEM_TARGET        = $NWCHEM_TARGET"
echo "ARMCI_NETWORK        = $ARMCI_NETWORK"
echo "USE_MPI              = $USE_MPI"
echo "LIBMPI               = $LIBMPI"
echo "MPI_LOC              = $MPI_LOC"
echo "MPI_INCLUDE          = $MPI_INCLUDE"
echo "MPI_LIB              = $MPI_LIB"
echo "LD_LIBRARY_PATH      = $LD_LIBRARY_PATH"
echo "NWCHEM_MODULES       = $NWCHEM_MODULES"
echo "USE_NOFSCHECK        = $USE_NOFSCHECK"
echo "USE_NOIO             = $USE_NOIO"
echo "LIB_DEFINES          = $LIB_DEFINES"
echo "MRCC_METHODS         = $MRCC_METHODS"
echo "PYTHONHOME           = $PYTHONHOME"
echo "PYTHONVERSION        = $PYTHONVERSION"
echo "PYTHON_EXE           = $PYTHON_EXE"
echo "USE_PYTHONCONFIG     = $USE_PYTHONCONFIG"
echo "USE_PYTHON64         = $USE_PYTHON64"
echo "PYTHONLIBTYPE        = $PYTHONLIBTYPE"
echo "PYTHONCONFIGDIR      = $PYTHONCONFIGDIR"
echo "LARGE_FILES          = $LARGE_FILES"
echo "USE_NOFSCHECK        = $USE_NOFSCHECK"
echo "NWCHEM_EXECUTABLE    = $NWCHEM_EXECUTABLE"
echo "NWCHEM_BASIS_LIBRARY = $NWCHEM_BASIS_LIBRARY"
echo "NWCHEM_NWPW_LIBRARY  = $NWCHEM_NWPW_LIBRARY"
echo "BLASOPT              = $BLASOPT"
