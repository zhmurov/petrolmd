 #!/bin/bash

# Set variables
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

if [[ $forcefield == "charmm36" ]]; then
    water=tip3p
elif [[ $forcefield == "trappeua" ]]; then
    water=tip4p
else
    echo "Only charmm36 and trappeua forefields are supported by this script"
    set -e
fi

for name in $hydrocarbons; do

    echo "================================================================================"
    echo "  Working on ${name}"
    echo "================================================================================"

    # Create folder for the molecular system
    mkdir ${name}
    cd ${name}

    # Copy and preapre packmol input
    cp ${PETROLMD}/SurfaceTension/files/template_water.inp packmol.inp
    sed -i "s/NEWMOLECULENAME/${name}/g" packmol.inp
    sed -i "s/FORCEFIELD/${forcefield}/g" packmol.inp
    sed -i "s/WATER/${water}/g" packmol.inp

    # Copy and prepare topology file
    cp ${PETROLMD}/SurfaceTension/files/template_water.top topol.top
    sed -i "s/NEWMOLECULENAME/${name}/g" topol.top
    sed -i "s/FORCEFIELD/${forcefield}/g" topol.top
    sed -i "s/WATER/${water}/g" topol.top

    # Get the configuration files
    cp ${PETROLMD}/files/mdp-${forcefield}/*.mdp .

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
    $GMX energy -f nvt.edr -xvg none -b 5000 -e 10000 <<< $'Pres-XX\nPres-YY\nPres-ZZ\n#Surf*SurfTen\n\n' -o nvt.pressure.xvg > nvt.pressure.out

    # Compute the surface tension
    python3 ${PETROLMD}/SurfaceTension/ComputeSurfaceTension.py nvt.pressure.xvg > gamma.txt

    cd ..
done
