 #!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
PACKMOL=~/git/external/packmol/packmol

SYSTEM_NAME=octane_solv2
PETROLMD=~/git/artemzhmurov/petrolmd
FFHOME=~/git/artemzhmurov/charmm36

cp -r ${FFHOME}/toppar/ .
cp -r ${FFHOME}/coord/ .

cp ${PETROLMD}/Octane/files/*.inp .
cp ${PETROLMD}/Octane/files/*.top .
cp ${PETROLMD}/files/mdp-charmm36/* .

$PACKMOL < ${SYSTEM_NAME}.inp
 
$GMX editconf -f ${SYSTEM_NAME}.pdb -o ${SYSTEM_NAME}.gro -box 4 4 4
$GMX grompp -f em.mdp -c ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -o em.tpr
$GMX mdrun -deffnm em
$GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
$GMX mdrun -deffnm nvt
$GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
$GMX mdrun -deffnm npt -update gpu
$GMX grompp -f md_iso.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_iso.tpr
$GMX mdrun -deffnm md_iso -update gpu
