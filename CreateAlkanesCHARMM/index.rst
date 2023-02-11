Linear hydrocarbons (alkanes) 
=============================

Coordinates
-----------

Topology
--------

In case with benzene molecule we were lucky since there was a residue topology that we used to create a molecular topology for the molecule.
Alkanes are a chain of carbon atoms saturated with hydrogens.
All carbons are the same apart from the first one.


To get an idea on what atom types and parameters to use, let us look at ``ethers.rtp`` file.
The very first residue there is ``BUTA`` (butane), which should be a good template for the rest of the alkanes.

    .. code-block:: text

        [ BUTA ]
        ; butane
        [ atoms ]
            H11    HCA3A  0.0900   1
            H12    HCA3A  0.0900   1
            H13    HCA3A  0.0900   1
            C1    CC33A -0.2700   1
            H21    HCA2A  0.0900   2
            H22    HCA2A  0.0900   2
            C2    CC32A -0.1800   2
            H31    HCA2A  0.0900   3
            H32    HCA2A  0.0900   3
            C3    CC32A -0.1800   3
            H41    HCA3A  0.0900   4
            H42    HCA3A  0.0900   4
            H43    HCA3A  0.0900   4
            C4    CC33A -0.2700   4
        [ bonds ]
            H11    C1
            H12    C1
            H13    C1
            C1    C2
            H21    C2
            H22    C2
            C2    C3
            H31    C3
            H32    C3
            C3    C4
            H41    C4
            H42    C4
            H43    C4


    .. code-block:: text

        ; Methylene group in the middle of an alkane
        [ CH2 ]
        [ atoms ]
            C      CC32A -0.1800   1
            H1     HCA2A  0.0900   1
            H2     HCA2A  0.0900   1
        [ bonds ]
        C     H1
        C     H2
        -C     C
        C    +C

    .. code-block:: text

        ; Terminal methyl group for an alkane
        [ CH3 ]
        [ atoms ]
            C      CC33A -0.2700   1
            H1     HCA3A  0.0900   1
            H2     HCA3A  0.0900   1
            H3     HCA3A  0.0900   1
        [ bonds ]
        C     H1
        C     H2
        C     H3
        -C     C


Create topologies and coordinates for the alkanes in CHARMM36
-------------------------------------------------------------

Make sure to use patched CHARMM36 forcefield from here: https://gitlab.com/artemzhmurov/charmm36 . 

The following will create initial coordinates for alkanes, build the topologies and minimize coordinates using GROMACS. Make sure that you change the path to GROMACS and to this repo in the script before running.

    .. code-block:: shell

        bash ${PETROLMD}/CreateAlkanesCHARMM/create_topologies.sh

This will produce two folders: coord, with minimized coordinates in gro and pdb formats and toppar with itp files for the molecules.

Setting variables
-----------------

    .. code-block:: shell

        GMX=/usr/local/gromacs/bin/gmx
        PACKMOL=~/git/external/packmol/packmol
        PETROLMD=~/git/artemzhmurov/petrolmd/
        CHARMM36_HOME=~/git/artemzhmurov/charmm36
        Lx=20
        Ly=20
        Lz=20
        SYSTEM_NAME=yamburg_recomb

Running simulations
-------------------

    .. code-block:: shell

        cp -r ${CHARMM36_HOME}/toppar/ .
        cp -r ${CHARMM36_HOME}/coord/ .
        cp ${PETROLMD}/files/mdp-charmm36/*.mdp .
        ${PETROLMD}/build/CountNumMolecules/count_mols ${PETROLMD}/CountNumMolecules/files/atomic_weights.dat ${PETROLMD}/CountNumMolecules/files/${SYSTEM_NAME}.dat ${SYSTEM_NAME} ${Lx} ${Ly} ${Lz}
        $PACKMOL < ${SYSTEM_NAME}_packmol.inp
        $GMX editconf -f ${SYSTEM_NAME}.pdb -o ${SYSTEM_NAME}_box.gro -box 30 30 30 -noc
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}_box.gro -p ${SYSTEM_NAME}.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
        $GMX mdrun -deffnm npt -update gpu
        $GMX grompp -f md_iso.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_iso.tpr
        $GMX mdrun -deffnm md_iso -update gpu
