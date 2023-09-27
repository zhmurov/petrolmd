#!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
SYSTEM_NAME=C6H6
PETROLMD=~/git/artemzhmurov/petrolmd
FORCEFIELD_HOME=~/git/external/trappeua

mkdir tmp
cd tmp

mkdir hydrocarbons

cp -r ${FORCEFIELD_HOME}/trappeua.ff .

# Create coordinates file for the C6H6 molecule
${PETROLMD}/build/BenzeneTraPPE-UA/create_benzene_trappe-ua


name=${SYSTEM_NAME}

# Create topology and minimize the structure
cp ${PETROLMD}/files/em_vac.mdp em.mdp

# Make a topology for the molecule
$GMX pdb2gmx -f ${name}.gro -o ${name}.gro -p ${name}.top -ff trappeua -water tip4p

# Create a copy of the topology that can be included
cp ${name}.top ${name}.itp
# Remove the header
sed -i -n '/\[ moleculetype \]/,$p' ${name}.itp
# Remove pairs section (needed by TraPPE forcefield)
if grep -Fxq "[ pairs ]" ${name}.itp
then
    sed -i -n '1,/pairs/p;/angles/,$p' ${name}.itp
    sed -i '\[ pairs \]/d' ${name}.itp
else
    echo "Skipping ${name}.itp"
fi
# Remove dihedrals section (in TraPPE forcefield aromatic rings are rigid)
if grep -Fxq "[ dihedrals ]" ${name}.itp
then
    sed -i -n '1,/dihedrals/p;/; Include Position restraint file/,$p' ${name}.itp
    sed -i '\[ dihedrals \]/d' ${name}.itp
else
    echo "Skipping ${name}.itp"
fi
# Remove the footer
sed -i '/; Include Position restraint file/,$d' ${name}.itp
# Rename the molecule
sed -i "s/Other/${name}/g" ${name}.itp
# Combine topologies into one file
# cat ${name}.itp >> alkanes.itp
# Copy the topolgy to separate folder
cp ${name}.itp hydrocarbons/${name}.itp

# Create topology that include itp file we just created
cp ${name}.top ${name}_bck.top
cp ${PETROLMD}/BenzeneTraPPE-UA/files/template.top ${name}.top
sed -i "s/NEWMOLECULENAME/${name}/g" ${name}.top

$GMX editconf -f ${name}.gro -o ${name}.gro -d 0.1
$GMX editconf -f ${name}.gro -o ${name}.gro -box 100 100 100 -noc
$GMX grompp -f em.mdp -c ${name}.gro -p ${name}.top -o ${name}_em.tpr
$GMX mdrun -deffnm ${name}_em

cp ${name}_em.gro hydrocarbons/${name}.gro
$GMX editconf -f ${name}_em.gro -o ${name}_em.pdb
cp ${name}_em.pdb hydrocarbons/${name}.pdb

