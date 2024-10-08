#!/bin/bash
GMX=/usr/local/gromacs/bin/gmx
PETROLMD=~/git/artemzhmurov/petrolmd

# Prepare output folder
mkdir hydrocarbons

# Copy the energy minimization script
cp ${PETROLMD}/files/em_vac.mdp em.mdp
# Generate PDBs for alkanes
${PETROLMD}/build/CreateAlkanesTraPPE-UA/create_alkanes
# Copy oter PDBs
cp ${PETROLMD}/CreateAlkanesTraPPE-UA/files/PDBs/* .

for filename in *.pdb; do
    name=$(basename $filename .pdb)

    # Create topology

    # Run GROMACS to create top files for all the PDBs in the folder
    $GMX pdb2gmx -f $filename -o ${name}.gro -p ${name}.top -i ${name}_posre.itp -ff trappeua -water tip4p
    # Create a copy of the topology that can be included
    cp ${name}.top ${name}.itp
    # Remove the header
    sed -i -n '/\[ moleculetype \]/,$p' ${name}.itp
    # Remove pairs section (needed by TraPPE forcefield)
    if grep -Fxq "[ pairs ]" ${name}.itp
    then
        sed -i -n '1,/pairs/p;/angles/,$p' ${name}.itp
        sed -i '\[ pairs \]/d' ${name}.itp
        sed -i -n '1,/pairs/p;/angles/,$p' ${name}.top
        sed -i '\[ pairs \]/d' ${name}.top
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
    
    # Minimize energy
    
    # Move molecule so that all coordinates are positive
    $GMX editconf -f ${name}.pdb -o ${name}.gro -d 0.1
    # Create lare simulation box (~no PBC simulations)
    $GMX editconf -f ${name}.gro -o ${name}.gro -box 100 100 100 -noc
    # Configure and run GROMACS
    $GMX grompp -f em.mdp -c ${name}.gro -p ${name}.top -o ${name}_em.tpr
    $GMX mdrun -deffnm ${name}_em
    # Rename the system in the resulting file
    sed -i "s/Protein/${name}/g" ${name}_em.gro
    # Convert GRO to PDB
    $GMX editconf -f ${name}_em.gro -o ${name}_em.pdb
    # Copy the resulting coordinates
    cp ${name}_em.gro hydrocarbons/${name}.gro
    cp ${name}.pdb hydrocarbons/${name}.pdb
done