#!/bin/bash
./SingleGrain
if [ $? == 0 ] #checks if the command executed before was executed successfully.
then
echo 10 > Results.sim 
awk 'NR==3' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==4' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==5' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==6' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==7' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==8' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==9' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==10' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==11' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
awk 'NR==12' Slope.dat | awk -F , '{ print "TimeStep1 " $5 " 0.0002"}'>>Results.sim
else 
echo "running benchmark SingleGrain failed"
fi

