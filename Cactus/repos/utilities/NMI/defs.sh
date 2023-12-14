#! /bin/bash

# Global settings (for all stages of the tests)

# This file is a shell script that is sourced.

echo '[BEGIN ENV]'
env
echo '[END ENV]'



case "$NMI_component" in
    
    "WaveToy")
        # WaveToy
        
        # WaveToy does not use an option list
        OPTIONLIST=
        
        # The thorn list is supposed to be checked out by one of the
        # input scripts listed in the submission file
        THORNLIST=Utilities/trunk/NMI/WaveToy.th
        
        # Name of the Cactus executable
        CONFIGURATION=wavetoy
        ;;
    
    "EinsteinToolkit")
        # CIGR
        
        # The option list is supposed to be checked out by one of the
        # input scripts listed in the submission file
        OPTIONLIST=Utilities/trunk/NMI/einsteintoolkit.cfg
        
        # The thorn list is supposed to be checked out by one of the
        # input scripts listed in the submission file
        THORNLIST=manifest/trunk/einsteintoolkit.th
        
        # Name of the Cactus executable
        CONFIGURATION=einsteintoolkit
        ;;
    
    *)
        echo "NMI_component not set"
        exit 1
        ;;
esac
