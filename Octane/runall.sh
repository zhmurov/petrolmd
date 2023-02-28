 #!/bin/bash

## Set variables
GMX=/usr/local/gromacs/bin/gmx
PACKMOL=~/git/external/packmol/packmol

PETROLMD=~/git/artemzhmurov/petrolmd
CHARMM36FFHOME=~/git/artemzhmurov/charmm36
TRAPPEUAFFHOME=~/git/external/trappeua

# Get the force-fields
cp -r ${CHARMM36FFHOME}/charmm36.ff .
cp -r ${TRAPPEUAFFHOME}/trappeua.ff .

hydrocarbons="C4H10_ISO C4H10 C5H12_ISO C5H12 C6H14 C6H6 C7H16 C8H18 C9H20 C10H22 C11H24 C12H26 C13H28 C14H30 C15H32 C16H34 C17H36 C18H38 C19H40 C20H42"

forcefield=$1

for name in $hydrocarbons; do

    echo "================================================================================"
    echo "  Working on ${name}"
    echo "================================================================================"

    # Create folder for the molecular system
    mkdir ${name}
    cd ${name}

    # Copy and preapre packmol input
    cp ${PETROLMD}/Octane/files/template.inp packmol.inp
    sed -i "s/NEWMOLECULENAME/${name}/g" packmol.inp
    sed -i "s/FORCEFIELD/${forcefield}/g" packmol.inp

    # Copy and prepare topology file
    cp ${PETROLMD}/Octane/files/template.top topol.top
    sed -i "s/NEWMOLECULENAME/${name}/g" topol.top
    sed -i "s/FORCEFIELD/${forcefield}/g" topol.top

    # Get the configuration files
    cp ${PETROLMD}/files/mdp-${forcefield}/*.mdp .

    $PACKMOL < packmol.inp
    
    $GMX editconf -f conf.pdb -o conf.gro -box 20 20 20
    $GMX grompp -f em.mdp -c conf.gro -o em.tpr
    $GMX mdrun -deffnm em
    $GMX grompp -f nvt.mdp -c em.gro -o nvt.tpr
    $GMX mdrun -deffnm nvt
    $GMX grompp -f npt.mdp -c nvt.gro -o npt.tpr
    $GMX mdrun -deffnm npt -update gpu
    $GMX grompp -f md_iso.mdp -c npt.gro -o md_iso.tpr
    $GMX mdrun -deffnm md_iso -update gpu

    cd ..

done
