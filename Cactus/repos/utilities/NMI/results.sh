#! /bin/bash

set -e                          # Output commands
set -x                          # Abort on errors

source Utilities/trunk/NMI/defs.sh



# Tar up all local files into a tarball "results.tar.gz", which will
# then be copied back to the submission host, where it is available
# for download over the web.  This is mostly for debugging, and could
# (or should) be disabled in the future.
tar czf results.tar.gz .

echo "Done."
