 #!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
PACKMOL=~/git/external/packmol/packmol

SYSTEM_NAME=octane_solv
PETROLMD=~/git/artemzhmurov/petrolmd
FFHOME=~/git/artemzhmurov/charmm36

cp -r ${FFHOME}/charmm36.ff .

cp ${PETROLMD}/SurfaceTension/files/*.inp .
cp ${PETROLMD}/SurfaceTension/files/*.top .
cp ${PETROLMD}/files/mdp-charmm36/* .

$PACKMOL < ${SYSTEM_NAME}.inp
 
$GMX editconf -f ${SYSTEM_NAME}.pdb -o ${SYSTEM_NAME}.gro -box 10 10 10
$GMX grompp -f em.mdp -c ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -o em.tpr
$GMX mdrun -deffnm em
$GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
$GMX mdrun -deffnm nvt
$GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
$GMX mdrun -deffnm npt -update gpu
$GMX grompp -f md_iso.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_iso.tpr
$GMX mdrun -deffnm md_iso -update gpu

$GMX energy -f md_iso.edr -xvg none -b 5000 <<< $'Pres-XX\nPres-YY\nPres-ZZ\n\n' -o md_iso.pressure.xvg > md_iso.pressure.out

