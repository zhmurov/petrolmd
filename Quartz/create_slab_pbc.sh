#!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
NX=10
NY=10
SYSTEM_NAME=slab_${NX}x${NY}
FFHOME=~/git/artemzhmurov/charmm36
PETROLMD=~/git/artemzhmurov/petrolmd

# Create coordinates file for the slab
${PETROLMD}/build/Quartz/create_quartz ${PETROLMD}/Quartz/files/input.xyz ${PETROLMD}/Quartz/files/crystal.dat ${SYSTEM_NAME}.gro yes ${NX} ${NY} 100.0

# Make a topology/coordinates pair for the entire molecule
cp ${SYSTEM_NAME}.top ${SYSTEM_NAME}_PBC.itp
sed -i -n '/\[ moleculetype \]/,$p' ${SYSTEM_NAME}_PBC.itp
sed -i '/; Include Position restraint file/,$d' ${SYSTEM_NAME}_PBC.itp
sed -i "s/Other/${SYSTEM_NAME}/g" ${SYSTEM_NAME}_PBC.itp
mkdir toppar
cp ${SYSTEM_NAME}.itp toppar/${SYSTEM_NAME}_PBC.itp
mkdir coord
cp ${SYSTEM_NAME}_em.gro coord/${SYSTEM_NAME}_PBC.gro
cp ${PETROLMD}/Quartz/files/sislab.top ${SYSTEM_NAME}_PBC.top
sed -i "s/NEWMOLECULENAME/${SYSTEM_NAME}_PBC/g" ${SYSTEM_NAME}_PBC.top

# Run test simulations
cp ${FFHOME}/specbond.dat .
$GMX pdb2gmx -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -ff charmm36 -water tip3p
$GMX solvate -cp ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top
cp ${PETROLMD}/files/mdp-charmm36/*.mdp .
$GMX grompp -f em.mdp -c ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top -o em.tpr
$GMX mdrun -deffnm em
$GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
$GMX mdrun -deffnm nvt
$GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
$GMX mdrun -deffnm npt -update gpu
$GMX grompp -f md_anis.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_anis.tpr
$GMX mdrun -deffnm md_anis -update gpu -nsteps 500000