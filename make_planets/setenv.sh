# set MESA_DIR to be the directory to which you downloaded MESA
# The directory shown is only an example and must be modified for your particular system.
export MESA_DIR=/home/daniele/codes/MESA/version_24.08.01/mesa-24.08.1

# set OMP_NUM_THREADS to be the number of cores on your machine
export OMP_NUM_THREADS=4

# you should have done this when you set up the MESA SDK
# The directory shown is only an example and must be modified for your particular system.
export MESASDK_ROOT=/home/daniele/codes/MESA/version_24.08.01/mesasdk
source $MESASDK_ROOT/bin/mesasdk_init.sh

# add shmesa (the MESA command line tool) to your PATH
export PATH=$PATH:$MESA_DIR/scripts/shmesa
