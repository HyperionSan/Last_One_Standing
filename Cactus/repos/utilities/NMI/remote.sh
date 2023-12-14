#! /bin/bash

set -e                          # Output commands
set -x                          # Abort on errors

source Utilities/trunk/NMI/defs.sh

# Should we use gmake or make?
MAKE=$(gmake --help > /dev/null 2>&1 && echo gmake || echo make)



case "$_NMI_TASKNAME" in
    
    "doc") 
        # TODO: Latex is not installed; we cannot build the
        # documentation on the build ndoes.  The documentation is
        # instead built only once, on the fetch node, after all other
        # tests have finished.
        
        # # Build documentation
        # pushd Cactus
        # $MAKE AllDoc
        ;;
    
    "configure")
        # Configure
        pushd Cactus
        $MAKE $CONFIGURATION-config ${OPTIONLIST:+options=../$OPTIONLIST} THORNLIST=../$THORNLIST <<EOF
yes
EOF
        ;;

    "build")
        # Build
        pushd Cactus
        $MAKE $CONFIGURATION -j2
        ;;

    "test")
        # Test
        pushd Cactus
        $MAKE $CONFIGURATION-testsuite <<EOF



EOF
        ;;
    
esac

echo "Done."
