#!/bin/bash
#------------------------------------------------------------------------------#
# CTS-1 Environment setups
#------------------------------------------------------------------------------#

export VENDOR_DIR=/usr/projects/draco/vendors

# symlinks will be generated for each machine that point to the correct
# installation directory.
if [[ `df | grep yellow | grep -c jayenne` -gt 0 ]]; then
  export siblings="snow badger kodiak grizzly"
else
  export siblings="fire ice cyclone"
fi

# The following toolchains will be used when releasing code
export environments="intel1802env intel1704env gcc740env"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

# SLURM
avail_queues=`sacctmgr -np list assoc user=$LOGNAME | sed -e 's/.*|\(.*dev.*\|.*access.*\)|.*/\1/' | sed -e 's/|.*//'`
case $avail_queues in
  *access*) access_queue="-A access --qos=access" ;;
  *dev*) access_queue="--qos=dev" ;;
esac
export access_queue

# Special setup for CTS-1: replace the 'latest' symlink
(cd /usr/projects/$package; if [[ -L latest ]]; then rm latest; fi; ln -s $source_prefix latest)

#------------------------------------------------------------------------------#
# Specify environments (modules)
#------------------------------------------------------------------------------#

if ! [[ $ddir ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (cts1-env.sh)"
  exit 1
fi

case $ddir in

  #------------------------------------------------------------------------------#
  draco-7_2* | draco-7_3*)
    function intel1802env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/18.0.2 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.4.0 gsl"
      run "module load mkl metis ndi csk qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    function intel1704env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/17.0.4 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.4.0 gsl"
      run "module load mkl metis ndi csk qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    function gcc740env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load gcc/7.4.0 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.4.0 gsl"
      run "module load mkl metis ndi csk qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    ;;


  #------------------------------------------------------------------------------#
  draco-6_25_0 | draco-7_0_0 | draco-7_1_0)
    function intel1802env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/18.0.2 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.3.0 gsl"
      run "module load mkl metis ndi csk qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    function intel1704env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/17.0.4 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.3.0 gsl"
      run "module load mkl metis ndi csk qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    function gcc640env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load gcc/6.4.0 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.3.0 gsl"
      run "module load mkl metis ndi qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    ;;

#------------------------------------------------------------------------------#
  draco-6_23_0 )

    function intel1704env()
    {
      run "module purge"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff"
      run "module load intel/17.0.4 openmpi/2.1.2"
      run "module load random123 eospac/6.2.4 gsl"
      run "module load mkl metis ndi csk"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }

    function intel1701env()
    {
      run "module purge"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff"
      run "module load intel/17.0.1 openmpi/1.10.5"
      run "module load random123 eospac/6.2.4 gsl"
      run "module load mkl metis ndi csk"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }

    function gcc640env()
    {
      run "module purge"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff"
      run "module load gcc/6.4.0 openmpi/2.1.2"
      run "module load random123 eospac/6.2.4 gsl"
      run "module load mkl metis ndi"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    ;;

esac

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##