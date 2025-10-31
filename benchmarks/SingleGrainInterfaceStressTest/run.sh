#!/bin/bash
./SingleGrainStressTest
if [ $? == 0 ]
then
    echo "Running benchmark interface stress test successful"
else 
    echo "Running benchmark interface stress test failed"
fi

