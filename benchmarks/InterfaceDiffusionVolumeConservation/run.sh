h!/bin/bash
./InterfaceDiffusionVolumeConservation
if [ $? == 0 ]
then
    echo "Running benchmark interface diffusion volume conservation successful"
else 
    echo "Running benchmark interface diffusion volume conservation failed"
fi

