#!/bin/bash
./InterfaceDiffusionJunction
if [ $? == 0 ]
then
    echo "Running benchmark interface diffusion junction successful"
else
    echo "Running benchmark interface diffusion junction failed"
fi

