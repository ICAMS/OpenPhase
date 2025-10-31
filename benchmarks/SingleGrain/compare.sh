#!/bin/bash
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color
paste Results.ref   Results.sim > compare.help #
awk '{n++; if(n>1){s=$1; x1=$2;  x2=$5; tol=$3; diff=x2-x1; if(diff<0.0)
diff=-diff; if (diff == 0) print s,x1,x2,tol, "EXACT"; else if (diff<tol) print s,x1,x2,tol, "INEXACT"; else print s,
x1,x2,tol, "WRONG"; }}' compare.help >compare.log
grep WRONG compare.log >/dev/null
if [ $? -eq 1 ] ; then
        echo -e "RESULTS: ${GREEN}OK${NC}" 
        echo -e "RESULTS: OK" >> compare.log
        grep INEXACT compare.log >/dev/null
        if [ $? -eq 0 ] ; then
            echo -e "${YELLOW}WARNING: ${NC}Results are within tolerance but have changed from a previous run:" 
            echo -e "WARNING: Results are within tolerance but have changed from a previous run:" >> compare.log
            grep INEXACT compare.log
        fi
        else
        echo -e "RESULTS: ${RED}FAILED${NC}" 
        echo -e "RESULTS: FAILED" >> compare.log
        grep WRONG compare.log
        fi
rm -f compare.help
