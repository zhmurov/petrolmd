#!/bin/bash
GMX=/usr/local/gromacs/bin/gmx

mkdir coord
cp ${PETROLMD}/files/em.mdp .

#for filename in *.pdb; do
filename=C1H4.pdb
    name=$(basename $filename .pdb)
    $GMX pdb2gmx -f ${name}.gro -o ${name}.gro -p ${name}.top -i ${name}_posre.itp -ff charmm36-jul2021 -water tip3p
    $GMX grompp -f em.mdp -c ${name}.gro -p ${name}.top -o ${name}_em.tpr
    $GMX mdrun -deffnm ${name}_em
    cp ${name}_em.gro coord/${name}.gro
#done