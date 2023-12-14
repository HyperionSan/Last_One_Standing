#! /bin/bash

set -e                          # Output commands
set -x                          # Abort on errors

source Utilities/trunk/NMI/defs.sh



# Copy the task list from its subdirectory to the current directory
cp Utilities/trunk/NMI/tasklist.nmi .



# Check out Cactus
Utilities/trunk/Scripts/GetComponents --anonymous --verbose $THORNLIST

echo "Done."
