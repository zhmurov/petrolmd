 #!/bin/bash

# Set variables
GMX=/usr/local/gromacs/bin/gmx
PACKMOL=~/git/external/packmol/packmol

PETROLMD=~/git/artemzhmurov/petrolmd
CHARMM36FFHOME=~/git/artemzhmurov/charmm36
TRAPPEUAFFHOME=~/git/external/trappeua

NUMMOLECULES=10000

LX=20
LY=20
LZ=12

let "LXA=${LX}*10"
let "LYA=${LY}*10"
let "LZA=${LZ}*10"

LZ2=100

let "LZ2A=${LZ2}*10"

# Get the force-fields
cp -r ${CHARMM36FFHOME}/charmm36.ff .
cp -r ${TRAPPEUAFFHOME}/trappeua.ff .

#hydrocarbons="C4H10_ISO C4H10 C5H12_ISO C5H12 C6H14 C6H6 C7H16 C8H18 C9H20 C10H22 C11H24 C12H26 C13H28 C14H30 C15H32 C16H34 C17H36 C18H38 C19H40 C20H42"
hydrocarbons="C10H22"

forcefield=$1

for name in $hydrocarbons; do

    echo "================================================================================"
    echo "  Working on ${name}"
    echo "================================================================================"

    # Create folder for the molecular system
    mkdir ${name}
    cd ${name}

    # Copy and preapre packmol input
    cp ${PETROLMD}/SurfaceTension/files/template_gas.inp packmol.inp
    sed -i "s/NEWMOLECULENAME/${name}/g" packmol.inp
    sed -i "s/FORCEFIELD/${forcefield}/g" packmol.inp
    sed -i "s/NUMMOLECULES/${NUMMOLECULES}/g" packmol.inp
    sed -i "s/LX/${LXA}/g" packmol.inp
    sed -i "s/LY/${LYA}/g" packmol.inp
    sed -i "s/LZ2/${LZ2A}/g" packmol.inp
    sed -i "s/LZ/${LZA}/g" packmol.inp

    # Copy and prepare topology file
    cp ${PETROLMD}/SurfaceTension/files/template_gas.top topol.top
    sed -i "s/NEWMOLECULENAME/${name}/g" topol.top
    sed -i "s/FORCEFIELD/${forcefield}/g" topol.top
    sed -i "s/WATER/${water}/g" topol.top
    sed -i "s/NUMMOLECULES/${NUMMOLECULES}/g" topol.top

    # Get the configuration files
    cp ${PETROLMD}/files/mdp-${forcefield}/*.mdp .

    # Create initial structure with Packmol
    $PACKMOL < packmol.inp
    
    # Configure, energy-minimize, equilibrate in GROMACS
    $GMX editconf -f conf.pdb -o conf.gro -box ${LX} ${LY} ${LZ2}
    $GMX grompp -f em.mdp -c conf.gro -o em.tpr
    $GMX mdrun -deffnm em
    $GMX grompp -f nvt.mdp -c em.gro -o nvt.tpr
    $GMX mdrun -deffnm nvt
    $GMX grompp -f npt.mdp -c nvt.gro -o npt.tpr
    $GMX mdrun -deffnm npt -update gpu -nsteps 500000

    # Production run    
    $GMX grompp -f md_nvt.mdp -c npt.gro -o md_nvt.tpr
    $GMX mdrun -deffnm md_nvt -update gpu

    # Get the pressure tensor components
    $GMX energy -f md_nvt.edr -xvg none -b 5000 -e 10000 <<< $'Pres-XX\nPres-YY\nPres-ZZ\n#Surf*SurfTen\n\n' -o md_nvt.pressure.xvg > md_nvt.pressure.out

    # Compute the surface tension
    python3 ${PETROLMD}/SurfaceTension/ComputeSurfaceTension.py md_nvt.pressure.xvg > gamma.txt

    cd ..
done
