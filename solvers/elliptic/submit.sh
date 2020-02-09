# LSF directives                                                                                                                                                               
#BSUB -P CSC262
#BSUB -J paranumal
#BSUB -o bp5.o%J
#BSUB -W 1:00

mesh_file=$1
N=$2
maxiter=$3
nodes=$4
ppn=$5

cd /gpfs/alpine/nfi114/proj-shared/thilina/bp_paper_libp_runs/libparanumal/solvers/elliptic
export OCCA_CACHE_DIR=`pwd`/.occa

# 1 First run this script to generate  OCCA kernels for only once
#./BP  ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}

# 2.After OCCA kernels are built, following script can be  used
jsrun -n${nodes} -r1 -a${ppn} -c${ppn} -g${ppn} ./BP  ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter} ${ppn}
