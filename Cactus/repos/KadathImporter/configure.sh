#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors

################################################################################
# Search
################################################################################

#if [ -z "${NETCDF_DIR}" ]; then
#    echo "BEGIN MESSAGE"
#    echo "Boost selected, but BOOST_DIR not set. Checking some places..."
#    echo "END MESSAGE"
#
#    FILES="include/boost/array.hpp include/boost/preprocessor/repetition/repeat.hpp include/boost/math/special_functions/gamma.hpp"
#    DIRS="/usr /usr/local /usr/local/boost /usr/local/packages/boost /usr/local/apps/boost ${HOME} c:/packages/boost"
#    for dir in $DIRS; do
#        BOOST_DIR="$dir"
#        for file in $FILES; do
#            if [ ! -r "$dir/$file" ]; then
#                unset BOOST_DIR
#                break
#            fi
#        done
#        if [ -n "$BOOST_DIR" ]; then
#            break
#        fi
#    done
#
#    if [ -z "$BOOST_DIR" ]; then
#        echo "BEGIN MESSAGE"
#        echo "Boost not found"
#        echo "END MESSAGE"
#    else
#        echo "BEGIN MESSAGE"
#        echo "Found Boost in ${BOOST_DIR}"
#        echo "END MESSAGE"
#    fi
#fi



################################################################################
# Build
################################################################################

#if [ -z "${BOOST_DIR}"                                                  \
#     -o "$(echo "${BOOST_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
#then
#    echo "BEGIN MESSAGE"
#    echo "Building Boost..."
#    echo "END MESSAGE"
#
#    # Set locations
#    THORN=Boost
#    #NAME=boost_1_54_0
#    NAME=boost_1_55_0
#    SRCDIR=$(dirname $0)
#    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
#    if [ -z "${BOOST_INSTALL_DIR}" ]; then
#        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
#    else
#        echo "BEGIN MESSAGE"
#        echo "Installing Boost into ${BOOST_INSTALL_DIR}"
#        echo "END MESSAGE"
#        INSTALL_DIR=${BOOST_INSTALL_DIR}
#    fi
#    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
#    BOOST_DIR=${INSTALL_DIR}
#
#    if [ -e ${DONE_FILE} -a ${DONE_FILE} -nt ${SRCDIR}/dist/${NAME}.tar.gz \
#                         -a ${DONE_FILE} -nt ${SRCDIR}/configure.sh ]
#    then
#        echo "BEGIN MESSAGE"
#        echo "Boost has already been built; doing nothing"
#        echo "END MESSAGE"
#    else
#        echo "BEGIN MESSAGE"
#        echo "Building Boost"
#        echo "END MESSAGE"
#
#        # Build in a subshell
#        (
#        exec >&2                # Redirect stdout to stderr
#        if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
#            set -x              # Output commands
#        fi
#        set -e                  # Abort on errors
#        cd ${SCRATCH_BUILD}
#
#        # Set up environment
#        unset LIBS
#        if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
#            export OBJECT_MODE=64
#        fi
#
#        echo "Boost: Preparing directory structure..."
#        mkdir build external done 2> /dev/null || true
#        rm -rf ${BUILD_DIR} ${INSTALL_DIR}
#        mkdir ${BUILD_DIR} ${INSTALL_DIR}
#
#        echo "Boost: Unpacking archive..."
#        pushd ${BUILD_DIR}
#        ${TAR?} xzf ${SRCDIR}/dist/${NAME}.tar.gz
#
#        echo "Boost: Configuring..."
#        cd ${NAME}
#        # Could also build with MPI instead
#        #./bootstrap.sh --prefix=${BOOST_DIR} --without-mpi
#        # Could also make MPI dependency optional
#        echo 'using mpi ;' > user-config.jam
#        ./bootstrap.sh --prefix=${BOOST_DIR}
#
#        echo "Boost: Building..."
#        makeflags=$(echo ${MAKE} | grep -e '-j' | perl -p -e 's/^.* -j *([0-9]+).*$/-j\1/')
#        ./b2 ${makeflags} link=static || true
#        #./b2 ${makeflags} || true
#
#        echo "Boost: Installing..."
#        ./b2 ${makeflags} install link=static || true
#        #./b2 ${makeflags} install || true
#        popd
#
#        echo "Boost: Cleaning up..."
#        rm -rf ${BUILD_DIR}
#
#        date > ${DONE_FILE}
#        echo "Boost: Done."
#        )
#
#        if (( $? )); then
#            echo 'BEGIN ERROR'
#            echo 'Error while building Boost. Aborting.'
#            echo 'END ERROR'
#            exit 1
#        fi
#    fi
#
#fi



################################################################################
# Configure Cactus
################################################################################


# Set include and lib dir
#KADATH_DIR=/home/astro/most/scratch/stuff/kadath
KADATH_INC_DIRS="${KADATH_DIR}/include"
KADATH_INC_DIRS="${KADATH_INC_DIRS} ${KADATH_DIR}/include/Kadath_point_h"
KADATH_LIB_DIRS="${KADATH_DIR}/lib"

# add -lVc if KADATH_LIBS is not set
if [ -z "${KADATH_LIBS}" ]; then
    KADATH_LIBS="kadath"
fi

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "HAVE_KADATH     = 1"
echo "KADATH_DIR      = ${KADATH_DIR}"
echo "KADATH_INC_DIRS = ${KADATH_INC_DIRS}"
echo "KADATH_LIB_DIRS = ${KADATH_LIB_DIRS}"
echo "KADATH_LIBS     = ${KADATH_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(KADATH_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(KADATH_LIB_DIRS)'
echo 'LIBRARY           $(KADATH_LIBS)'
