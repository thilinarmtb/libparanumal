#!/bin/bash                                                                                                                                                                     

function log2() {
  local input=$1
  local output=0

  while ((input>1)); do
    ((input/=2, output +=1))
  done

  echo $output
}

maxiter=1000;
mesh_base=../../meshes/box;
mesh_extension=.msh;
starting_elements=5
ppn=6

for nodes in 1; do
for elem in 4; do
#for nodes in 4; do
#for elem in 1 2 3 4 5 6 7 8 9 10 11 12; do
for N in 12; do
  ((mesh_id=starting_elements+elem))
  meshfile=${mesh_base}${mesh_id}${mesh_extension}
  echo "Running ${meshfile} on ${nodes} node(s) at order ${N}"
  bsub -nnodes ${nodes} -W 1:00 -P nfi114 -J paranumal -o bp5.o%J \
        ./submit.sh ${meshfile} ${N} ${maxiter} ${nodes} ${ppn}
  sleep 3;
done;
done;
done;
