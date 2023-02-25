 #!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
PACKMOL=~/git/external/packmol/packmol

SYSTEM_NAME=C8H18
PETROLMD=~/git/artemzhmurov/petrolmd
FFHOME=~/git/artemzhmurov/charmm36

# Get the force-field
cp -r ${FFHOME}/charmm36.ff .

hydrocarbons="C4H10_ISO C4H10 C5H12_ISO C5H12 C6H14 C6H6 C7H16 C8H18 C9H20 C10H22 C11H24 C12H26 C13H28 C14H30 C15H32 C16H34 C17H36 C18H38 C19H40 C20H42"

for name in $hydrocarbons; do

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
    cp ${PETROLMD}/files/mdp-charmm36/*.mdp .

    # Create initial structure with Packmol
    $PACKMOL < packmol.inp
    
    # Configure, energy-minimize, equilibrate in GROMACS
    $GMX editconf -f conf.pdb -o conf.gro -box 5 5 40
    $GMX grompp -f em.mdp -c conf.gro -o em.tpr
    $GMX mdrun -deffnm em
    $GMX grompp -f nvt.mdp -c em.gro -o nvt.tpr
    $GMX mdrun -deffnm nvt
    $GMX grompp -f npt.mdp -c nvt.gro -o npt.tpr
    $GMX mdrun -deffnm npt -update gpu

    # Production run    
    $GMX grompp -f md_iso.mdp -c npt.gro -o md_iso.tpr
    $GMX mdrun -deffnm md_iso -update gpu

    # Get the pressure tensor components
    $GMX energy -f md_iso.edr -xvg none -b 1000 -e 2000 <<< $'Pres-XX\nPres-YY\nPres-ZZ\n\n' -o md_iso.pressure.xvg > md_iso.pressure.out

    cd ..
done
