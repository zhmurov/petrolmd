 #!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
PACKMOL=~/git/external/packmol/packmol

SYSTEM_NAME=C8H18
PETROLMD=~/git/artemzhmurov/petrolmd
FFHOME=~/git/artemzhmurov/charmm36

# Get the force-field
cp -r ${FFHOME}/charmm36.ff .

for name in C1H4 C2H6 C3H8 C4H10_ISO C4H10 C5H12_ISO C5H12 C6H14 C6H6 C7H16 C8H18 C9H20 C10H22 C11H24 C12H26 C13H28 C14H30 C15H32 C16H34 C17H36 C18H38 C19H40 C20H42; do

    # Create folder for the molecular system
    mkdir ${name}
    cd ${name}

    # Copy and preapre packmol input
    cp ${PETROLMD}/SurfaceTension/files/template.inp packmol.inp
    sed -i "s/NEWMOLECULENAME/${name}/g" packmol.inp

    # Copy and prepare topology file
    cp ${PETROLMD}/SurfaceTension/files/template.top topol.top
    sed -i "s/NEWMOLECULENAME/${name}/g" topol.top

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
    $GMX mdrun -deffnm nvt -nsteps 1000000

    # Get the pressure tensor components
    $GMX energy -f nvt.edr -xvg none -b 1000 -e 2000 <<< $'Pres-XX\nPres-YY\nPres-ZZ\n\n' -o nvt.pressure.xvg > nvt.pressure.out

    cd ..
done
