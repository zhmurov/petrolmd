#!/bin/bash
GMX=/usr/local/gromacs/bin/gmx

mkdir coord
#cp ${PETROLMD}/files/em.mdp .
cp ~/git/artemzhmurov/petrolmd/CreateAlkanesCHARMM/files/em.mdp .
~/git/artemzhmurov/petrolmd/build/CreateAlkanesCHARMM/create_alkanes_charmm

for filename in *.pdb; do
#filename=C8H18.pdb
    name=$(basename $filename .pdb)
    $GMX editconf -f ${name}.gro -o ${name}.gro -d 0.1
    $GMX editconf -f ${name}.gro -o ${name}.gro -box 100 100 100 -noc
    $GMX pdb2gmx -f ${name}.gro -o ${name}.gro -p ${name}.top -i ${name}_posre.itp -ff charmm36-jul2021 -water tip3p
    $GMX grompp -f em.mdp -c ${name}.gro -p ${name}.top -o ${name}_em.tpr
    $GMX mdrun -deffnm ${name}_em
    $GMX mdrun -deffnm ${name}_em -c ${name}.pdb
    cp ${name}_em.gro coord/${name}.gro
    cp ${name}.pdb coord/${name}.pdb
done