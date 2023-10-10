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

systems="C6H6"
# temperatures="260 270 280 290 300 310 320"
temperatures="260"

forcefield=$1

for name in $systems; do

    echo "================================================================================"
    echo "  Working on ${name}"
    echo "================================================================================"

    # Create folder for the molecular system
    mkdir ${name}
    cd ${name}
    cp -r ${TRAPPEUAFFHOME}/trappeua.ff .

    # Copy and preapre packmol input
    cp ${PETROLMD}/BenzeneTraPPE-UA/files/template.inp packmol.inp
    sed -i "s/NEWMOLECULENAME/${name}/g" packmol.inp
    sed -i "s/FORCEFIELD/${forcefield}/g" packmol.inp

    # Copy and prepare topology file
    cp ${PETROLMD}/BenzeneTraPPE-UA/files/template_1000.top topol.top
    sed -i "s/NEWMOLECULENAME/${name}/g" topol.top
    sed -i "s/FORCEFIELD/${forcefield}/g" topol.top

    # Get the configuration files
    cp ${PETROLMD}/BenzeneTraPPE-UA/files/em.mdp .

    $PACKMOL < packmol.inp
    
    $GMX editconf -f conf.pdb -o conf.gro -box 10 10 10
    $GMX grompp -f em.mdp -c conf.gro -o em.tpr -maxwarn 1
    $GMX mdrun -deffnm em

    for temperature in $temperatures; do

        # Get the configuration files and set the temperature there
        cp ${PETROLMD}/BenzeneTraPPE-UA/files/nvt_T.mdp nvt_${temperature}.mdp
        sed -i "s/TEMPERATURE/${temperature}/g" nvt_${temperature}.mdp
        cp ${PETROLMD}/BenzeneTraPPE-UA/files/npt_T.mdp npt_${temperature}.mdp
        sed -i "s/TEMPERATURE/${temperature}/g" npt_${temperature}.mdp
        cp ${PETROLMD}/BenzeneTraPPE-UA/files/md_iso_T.mdp md_iso_${temperature}.mdp
        sed -i "s/TEMPERATURE/${temperature}/g" md_iso_${temperature}.mdp
        
        $GMX grompp -f nvt_${temperature}.mdp -c em.gro -o nvt_${temperature}.tpr -maxwarn 1
        $GMX mdrun -deffnm nvt_${temperature}
        $GMX grompp -f npt_${temperature}.mdp -c nvt_${temperature}.gro -o npt_${temperature}.tpr -maxwarn 1
        $GMX mdrun -deffnm npt_${temperature} -update gpu

        $GMX grompp -f md_iso_${temperature}.mdp -c npt_${temperature}.gro -o md_iso_${temperature}.tpr -maxwarn 1
        $GMX mdrun -deffnm md_iso_${temperature} -update gpu

    done

    cd ..

done
