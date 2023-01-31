#!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
SYSTEM_NAME=slab
FFHOME=~/git/artemzhmurov/charmm36
PETROLMD=~/git/artemzhmurov/petrolmd

# Create coordinates file for the slab
${PETROLMD}/build/Quartz/create_quartz

# Create topology and minimize the structure
cp ${FFHOME}/specbond.dat .
cp ${PETROLMD}/files/em_vac.mdp em.mdp
$GMX pdb2gmx -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -ff charmm36 -water tip3p
$GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -d 0.1
$GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -box 100 100 100 -noc
$GMX grompp -f em.mdp -c ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -o ${SYSTEM_NAME}_em.tpr
$GMX mdrun -deffnm ${SYSTEM_NAME}_em


# Make a topology/coordinates pair for the entire molecule
cp ${SYSTEM_NAME}.top ${SYSTEM_NAME}.itp
sed -i -n '/\[ moleculetype \]/,$p' ${SYSTEM_NAME}.itp
sed -i '/; Include Position restraint file/,$d' ${SYSTEM_NAME}.itp
sed -i "s/Other/${SYSTEM_NAME}/g" ${SYSTEM_NAME}.itp
mkdir toppar
cp ${SYSTEM_NAME}.itp toppar/${SYSTEM_NAME}.itp
mkdir coord
cp ${SYSTEM_NAME}_em.gro coord/${SYSTEM_NAME}.gro
cp ${PETROLMD}/Quartz/files/sislab.top ${SYSTEM_NAME}.top
sed -i "s/NEWMOLECULENAME/${SYSTEM_NAME}/g" ${SYSTEM_NAME}.top

# Run test simulations
$GMX editconf -f coord/${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -box 10 10 10 -noc
$GMX solvate -cp ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top
cp ${PETROLMD}/files/mdp-charmm36/*.mdp .
$GMX grompp -f em.mdp -c ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top -o em.tpr
$GMX mdrun -deffnm em
$GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
$GMX mdrun -deffnm nvt
$GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
$GMX mdrun -deffnm npt -update gpu
$GMX grompp -f md_iso.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_iso.tpr
$GMX mdrun -deffnm md_iso -update gpu -nsteps 500000