 #!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
PACKMOL=~/git/external/packmol/packmol

SYSTEM_NAME=C8H18
PETROLMD=~/git/artemzhmurov/petrolmd
FFHOME=~/git/artemzhmurov/charmm36

# Get the force-field
cp -r ${FFHOME}/charmm36.ff .

# Create folder for the molecular system
mkdir ${SYSTEM_NAME}
cd ${SYSTEM_NAME}

# Copy and preapre packmol input
cp ${PETROLMD}/SurfaceTension/files/template.inp packmol.inp
sed -i "s/NEWMOLECULENAME/${SYSTEM_NAME}/g" packmol.inp

# Copy and prepare topology file
cp ${PETROLMD}/SurfaceTension/files/template.top topol.top
sed -i "s/NEWMOLECULENAME/${SYSTEM_NAME}/g" topol.top

# Get the configuration files
cp ${PETROLMD}/files/mdp-charmm36/em.mdp .
cp ${PETROLMD}/files/mdp-charmm36/nvt.mdp .

# Create initial structure with Packmol
$PACKMOL < packmol.inp
 
# Configure and energy-minimize GROMACS
$GMX editconf -f conf.pdb -o conf.gro -box 40 5 5
$GMX grompp -f em.mdp -c conf.gro -o em.tpr
$GMX mdrun -deffnm em

# Production run
$GMX grompp -f nvt.mdp -c em.gro -o nvt.tpr
$GMX mdrun -deffnm nvt -nsteps 5000000

# Get the pressure tensor components
$GMX energy -f nvt.edr -xvg none -b 5000 <<< $'Pres-XX\nPres-YY\nPres-ZZ\n\n' -o nvt.pressure.xvg > nvt.pressure.out

cd ..
