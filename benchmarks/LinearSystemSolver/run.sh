#!/bin/bash
./LinearSystemSolver
if [ $? == 0 ]
then
    echo "Running benchmark LinearSystemSolver successful"
else
    echo "Running benchmark LinearSystemSolver failed"
fi

