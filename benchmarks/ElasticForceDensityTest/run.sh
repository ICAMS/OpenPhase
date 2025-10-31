#!/bin/bash
./ElasticForceDensityTest
if [ $? == 0 ]
then
    echo "Running benchmark elastic force density successful"
else
    echo "Running benchmark elastic force density failed"
fi

