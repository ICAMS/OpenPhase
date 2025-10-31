#!/bin/bash
./LBCapillaryBridge
if [ $? == 0 ]
then
    echo "Running benchmark capillary bridge successful"
else
    echo "Running benchmark capillary bridge failed"
fi

