Linear hydrocarbons (alkanes) 
=============================

The aim of this section is to create coordinates and molecular topologies for all alkanes that may be found in the mixture.
Likely there are some already present in the force-field and these will be used as a source of parameters and atom types.

Coordinates
-----------

Alkanes are chains of carbon atoms saturated with hydrogens.

Let us first look at the butane topology file to extract all the necessary parameters for the geometry of alkanes.

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

Let us start with carbon atoms.
The first and last carbons are of type ``CC33A`` and have a partial charge -0.27.
This type is a carbon, that is covalently linked to one other carbon and three hydrogens.
The rest of carbons have type ``CC32A`` and partial charge of -0.18.
These are linked to two other carbons in the chain and have two hydrogens attached to them.

The hydrogens that are connected to the first and last carbon are of type ``HCA3A``, the rest are of type ``HCA2A``.
All hydrogens have partial charge of 0.9, which makes the total charge of the molecule zero.

Let us check the bonded parameters for these atom types in ``ffbonded.itp`` file.
We are interested in equilibrium bond distances and angles, since we are going to use these parameters when constructing the coordinates of the molecule.
Corresponding lines for bond distances in ``ffbonded.itp`` file are:

    .. code-block:: text

        [ bondtypes ]
        ...
        ;      i        j  func           b0           kb
        ...
        CC32A    HCA2A     1   0.11110000    258571.20 ; alkanes, 4/98
        CC33A    HCA3A     1   0.11110000    269449.60 ; alkanes, 4/98
        ...
        CC32A    CC32A     1   0.15300000    186188.00 ; alkanes, 3/92
        CC32A    CC33A     1   0.15280000    186188.00 ; alkanes, 3/92
        CC33A    CC33A     1   0.15300000    186188.00 ; alkanes, 3/92
        ...

Each line here contains parameters for one type of bond.
The later is defined by the atom types that participate in this bond, which are in the first two columns.
These are followed by the function type (1 for harmonic bond), equilibrium distance (``b0``) and spring constant (``kb``).
Most of the entries also have a comment which indicates where the parameters are taken from.

Here, we see that the equilibrium distances for C-H bonds are the same for two cases we have: 0.1111 nm.
There is a slight variation in C-C bond distances, with bonds at the terminals of the molecule slightly lower in length.
For simplicity, we are going to use the same value of 0.153 nm for all of these and will let the energy minimization algorithm to fix this small inconsistency.

The parameters for equilibrium bond angles can be taken from the following lines in ``ffbonded.itp`` file:

    .. code-block:: text

        [ angletypes ]
        ...
        ;      i        j        k  func       theta0       ktheta          ub0          kub

        ...
        HCA2A    CC32A    CC32A     5   110.100000   221.752000   0.21790000     18853.10 ; alkane, 4/98
        HCA2A    CC32A    CC33A     5   110.100000   289.532800   0.21790000     18853.10 ; alkane, 4/98
        ...
        HCA3A    CC33A    CC33A     5   110.100000   313.800000   0.21790000     18853.10 ; alkane, 4/98
        HCA2A    CC32A    HCA2A     5   109.000000   297.064000   0.18020000      4518.72 ; alkane, 3/92
        HCA3A    CC33A    HCA3A     5   108.400000   297.064000   0.18020000      4518.72 ; alkane, 3/92
        ...
        CC32A    CC32A    CC32A     5   113.600000   488.272800   0.25610000      9338.69 ; alkane, 3/92
        ...
        CC32A    CC32A    CC33A     5   115.000000   485.344000   0.25610000      6694.40 ; alkane, 3/92

To get the geometrical parameters of a chain

Topology
--------

In case with benzene molecule we were lucky since there was a residue topology that we used to create a molecular topology for the molecule.

All carbons are the same apart from the first one.


To get an idea on what atom types and parameters to use, let us look at ``ethers.rtp`` file.
The very first residue there is ``BUTA`` (butane), which should be a good template for the rest of the alkanes.




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
