#! /bin/bash

set -e                          # Output commands
set -x                          # Abort on errors

source Utilities/trunk/NMI/defs.sh

# Should we use gmake or make?
MAKE=$(gmake --help > /dev/null 2>&1 && echo gmake || echo make)



# Generate all documentation
pushd Cactus
$MAKE AllDoc

echo "Done."
