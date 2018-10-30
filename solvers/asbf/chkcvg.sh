#!/bin/bash

# Maximum polynomial degree.
DMAX=12

# Maximum no. of radial modes.
NMAX=12

BASE_SETUP="setups/setupQuad3D.rc"
TEMP_SETUP="setupcvgchk.rc"

> "${TEMP_SETUP}"

I=1
while [ ${I} -le $DMAX ] ; do
  J=1
  while [ ${J} -le $NMAX ] ; do
    sed "30 s/.*/${I}/g; 33 s/.*/${J}/g" < "${BASE_SETUP}" > "${TEMP_SETUP}"
    ERRORS=$(./asbfMain "${TEMP_SETUP}" | tail -n1)
    ERR_H1=$(echo "${ERRORS}" | cut -f 1 -d ' ')
    ERR_L2=$(echo "${ERRORS}" | cut -f 2 -d ' ')
    echo "${I},${J},${ERR_H1},${ERR_L2}"
    J=$((J + 1))
  done
  I=$((I + 1))
done

#rm -f "${TEMP_SETUP}"
