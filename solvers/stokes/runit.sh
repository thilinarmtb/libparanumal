#!/bin/bash

setupFile()
{
	cat <<EOF
[FORMAT]
1.0

##############################################
#################### MESH ####################
##############################################

[MESH FILE]
../../meshes/unitSquareQuad2D-5.msh

[DATA FILE]
data/stokesDirichletTestQuad2D.h

[MESH DIMENSION]
2

[ELEMENT TYPE]
4

[INTEGRATION TYPE]
GLL
#CUBATURE

################################################
#################### DEVICE ####################
################################################

[THREAD MODEL]
#OpenCL
CUDA
#Serial

[PLATFORM NUMBER]
0

[DEVICE NUMBER]
0

################################################
################ PROBLEM SETUP #################
################################################

[LAMBDA]
${1}

[VELOCITY POLYNOMIAL DEGREE]
5

[VELOCITY DISCRETIZATION]
CONTINUOUS

[PRESSURE POLYNOMIAL DEGREE]
4

[PRESSURE DISCRETIZATION]
CONTINUOUS

################################################
############ KRYLOV SOLVER OPTIONS #############
################################################

[KRYLOV SOLVER]
MINRES
#DQGMRES

[KRYLOV SOLVER ITERATION LIMIT]
1000

[KRYLOV SOLVER TOLERANCE]
1.0e-6

################################################
############ PRECONDITIONER OPTIONS ############
################################################

[PRECONDITIONER]
#NONE
#JACOBI
SCHURCOMPLEMENTBLOCKDIAG

[JACOBI BOOST PARAMETER]
1.0e-1

[VELOCITY BLOCK PRECONDITIONER]
#MULTIGRID
JACOBI

[PRESSURE BLOCK PRECONDITIONER]
#MASSMATRIX
MULTIGRID

################################################
########## VELOCITY MULTIGRID OPTIONS ##########
################################################

[VELOCITY MULTIGRID COARSENING]
HALFDOFS

[VELOCITY MULTIGRID SMOOTHER]
DAMPEDJACOBI+CHEBYSHEV

[VELOCITY MULTIGRID CHEBYSHEV DEGREE]
1

[VELOCITY MULTIGRID PARALMOND AGGREGATION STRATEGY]
DEFAULT

[VELOCITY MULTIGRID PARALMOND CYCLE]
KCYCLE

[VELOCITY MULTIGRID PARALMOND SMOOTHER]
DAMPEDJACOBI+CHEBYSHEV

[VELOCITY MULTIGRID PARALMOND CHEBYSHEV DEGREE]
1

################################################
########## PRESSURE MULTIGRID OPTIONS ##########
################################################

[PRESSURE MULTIGRID COARSENING]
HALFDOFS

[PRESSURE MULTIGRID SMOOTHER]
DAMPEDJACOBI+CHEBYSHEV

[PRESSURE MULTIGRID CHEBYSHEV DEGREE]
1

[PRESSURE MULTIGRID PARALMOND AGGREGATION STRATEGY]
DEFAULT

[PRESSURE MULTIGRID PARALMOND CYCLE]
KCYCLE

[PRESSURE MULTIGRID PARALMOND SMOOTHER]
DAMPEDJACOBI+CHEBYSHEV

[PRESSURE MULTIGRID PARALMOND CHEBYSHEV DEGREE]
1

################################################
################ OTHER OPTIONS #################
################################################

[VERBOSE]
TRUE
EOF
}

echo "LAMBDA      Iterations   Error"
echo "--------------------------------"
for LAMBDA in "0.0" "1.0e8" "1.0e7" "1.0e6" "1.0e5" "1.0e4" "1.0e3" "1.0e2" "1.0e1" "1.0e0" ; do
#for LAMBDA in "0.0" "1.0e8"  ; do
	setupFile "${LAMBDA}" > setup.rc
	./stokesMain setup.rc > output.log

	ITERS="$(grep iterations output.log | cut -d ' ' -f 4)"
	ERROR="$(printf "%.1e\n" "$(grep errxInf output.log | cut -d '=' -f 2)")"

	printf "%6.1e     %3d          %6.1e\n" "${LAMBDA}" "${ITERS}" "${ERROR}"
done
