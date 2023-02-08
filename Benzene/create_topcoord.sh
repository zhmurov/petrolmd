#!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
SYSTEM_NAME=C6H6
PETROLMD=~/git/artemzhmurov/petrolmd

# Create coordinates file for the C6H6 molecule
${PETROLMD}/build/Benzene/create_benzene

# Create topology and minimize the structure
cp ${PETROLMD}/files/em_vac.mdp em.mdp
$GMX pdb2gmx -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -ff charmm36 -water tip3p
$GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -d 0.1
$GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -box 100 100 100 -noc
$GMX grompp -f em.mdp -c ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -o ${SYSTEM_NAME}_em.tpr
$GMX mdrun -deffnm ${SYSTEM_NAME}_em
mkdir coord
cp ${SYSTEM_NAME}_em.gro coord/${SYSTEM_NAME}.gro
$GMX editconf -f ${SYSTEM_NAME}_em.gro -o ${SYSTEM_NAME}.pdb
cp ${SYSTEM_NAME}.pdb coord/${SYSTEM_NAME}.pdb

# Make a topology for the molecule
cp ${SYSTEM_NAME}.top ${SYSTEM_NAME}.itp
sed -i -n '/\[ moleculetype \]/,$p' ${SYSTEM_NAME}.itp
sed -i '/; Include Position restraint file/,$d' ${SYSTEM_NAME}.itp
sed -i "s/Other/${SYSTEM_NAME}/g" ${SYSTEM_NAME}.itp
mkdir toppar
cp ${SYSTEM_NAME}.itp toppar/${SYSTEM_NAME}.itp

# Example on how to use the molecule topology
cp ${PETROLMD}/Benzene/files/template.top ${SYSTEM_NAME}.top
sed -i "s/NEWMOLECULENAME/${SYSTEM_NAME}/g" ${SYSTEM_NAME}.top