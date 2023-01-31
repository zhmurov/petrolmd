Creating the coordinates and topology for a quartz slab of a given size
-----------------------------------------------------------------------

To construct coordinates for a slab of SiO2, we will need (1) coordinates for a single unit cell and (2) size of the unit cell in two dimensions. Note that these files should be in-sync with the topology of a single unit in the forcefield, i.e. the atom order should be the same.

    .. code-block:: text

        [ SiO ]
        ; 
        [ atoms ]
            SI1     SIH3  0.8000   0
            SI2       SI  1.0800   0
            SI3       SI  1.0800   0
            SI4       SI  1.0800   0
            O5     OSIE -0.5300   0
            O6     OSIE -0.5300   0
            O7     OSIE -0.5300   0
            O8     OSIE -0.5300   0
            O9     OSIE -0.5300   0
            O10     OSIE -0.5300   0
            O11     OSIE -0.5300   0
            O12     OSIE -0.5300   0
            O13     OSIE -0.5300   0
            O14     OSIE -0.5300   0
            O15     OSIE -0.5300   0
            O16     OSIH -0.5400   0
            O17     OSIE -0.5300   0
            O18     OSIE -0.5300   0
            O19     OSIE -0.5300   0
            O20     OSIE -0.5300   0
            SI21       SI  1.0800   0
            SI22       SI  1.0800   0
            SI23       SI  1.0800   0
            SI24       SI  1.0800   0
            SI25     SIH3  0.8000   0
            SI26       SI  1.0800   0
            SI27       SI  1.0800   0
            SI28       SI  1.0800   0
            O29     OSIE -0.5300   0
            O30     OSIE -0.5300   0
            O31     OSIE -0.5300   0
            O32     OSIH -0.5400   0
            O33     OSIE -0.5300   0
            O34     OSIE -0.5300   0
            O35     OSIE -0.5300   0
            O36     OSIE -0.5300   0
            H37     HSIA -0.1500   0
            H38     HSIO  0.3200   0
            H39     HSIA -0.1500   0
            H40     HSIO  0.3200   0
        [ bonds ]
            O5   SI1
            O9   SI1
            H37   SI1
            O6   SI2
            O10   SI2
            O13   SI2
            O18   SI2
            O7   SI3
            O11   SI3
            O14   SI3
            O19   SI3
            O15   SI4
            O20   SI4
            SI21    O5
            SI24    O8
            SI26   O10
            SI27   O11
            SI28   O12
            SI25   O13
            SI26   O14
            SI27   O15
            SI28   O16
            SI21   O17
            SI23   O19
            SI24   O20
            O29  SI21
            O33  SI21
            O30  SI22
            O34  SI22
            O31  SI23
            O35  SI23
            O32  SI24
            O36  SI24
            H39  SI25
            O29  SI26
            O34  SI26
            O36  SI28
            H40   O32
            O16   H38

To link units covalently into a solid slab, we are going to use the special bonds mechanisms in GROMACS. Another example when this mechanism is used is the disulphide bonds in proteins, or bonds in HEM in haemoglobin. The file in which the special bonds are listed is shared among all the forcefields and located in the root of the forcefield folder. As with the rest of such files, it can be overriden by local file with the same name (i.e. `specbonds.dat`), located in the simulation folder. In this file, the first line is the number of special bonds entries in the file, followed by the list of the bonds, one per line. The format of the lines is ``resA atomA nbondsA resB atomB nbondsB length newresA newresB``, where:

1. ``resA`` and ``resB`` are the names of two residues, connected by the bond.

2. ``atomA`` and ``AtomB`` are the names of the connected atoms.

3. ``nbondsA`` and ``nbondsB`` are the number of times each of the atoms can be connected via the special bonds.

4. ``length`` is a reference length of the bond. If the distance between two atoms ``A`` and ``B`` above is within 10% of this length, the bond between them will be added.

5. ``newresA`` and ``newresB`` is the name of the residues after they were connected by the special bond. Should be the same as ``resA`` and ``resB`` if the names are not to be changed.

The desc More on the special bonds file in GROMACS manual `pdb2gmx documentation <https://manual.gromacs.org/current/reference-manual/topologies/pdb2gmx-input-files.html#specbond>`_).

For the topology file below, the following special bonds are needed:

    .. code-block:: text

        SiO     SI1     1       SiO     O17     1       0.16    SiO     SiO
        SiO     SI4     1       SiO     O8      1       0.26    SiO     SiO
        SiO     SI4     1       SiO     O12     1       0.16    SiO     SiO
        SiO     SI22    1       SiO     O6      1       0.26    SiO     SiO
        SiO     SI22    1       SiO     O18     1       0.16    SiO     SiO
        SiO     SI23    1       SiO     O7      1       0.26    SiO     SiO
        SiO     SI25    1       SiO     O9      1       0.16    SiO     SiO
        SiO     SI25    1       SiO     O33     1       0.30    SiO     SiO
        SiO     SI27    1       SiO     O30     1       0.18    SiO     SiO
        SiO     SI27    1       SiO     O35     1       0.16    SiO     SiO
        SiO     SI28    1       SiO     O31     1       0.18    SiO     SiO

Note, that if these bonds are added on top of those listed in the topology file above, all the atoms are saturated with covalent bonds: all the silicone atoms have four bonds in total, all oxygens --- two. Also note, that somme of the target distances are quite far off from the forcefield target distance of 0.1698 nm (see the bond between atom types ``SI`` and ``OSIE`` in ``ffbonded.itp`` file of the forcefield). These distances are taken from the slab structure and used here to make sure that corresponding bonds will be added (i.e. will be within the 10% margin of the distance). Afterwards, the energy minimization algorithm should pull the respective atoms closer to one another or push them farther apart.

    .. code-block:: shell

        GMX=/usr/local/gromacs/bin/gmx
        SYSTEM_NAME=slab
        FFHOME=~/git/artemzhmurov/charmm36
        PETROLMD=~/git/artemzhmurov/petrolmd
        cp ${FFHOME}/specbond.dat .
        cp ${PETROLMD}/files/em_vac.mdp em.mdp
        ${PETROLMD}/build/Quartz/create_quartz
        $GMX pdb2gmx -f ${SYSTEM_NAME}.pdb -o ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -ff charmm36 -water tip3p
        $GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -d 0.1
        $GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -box 100 100 100 -noc
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -o ${SYSTEM_NAME}_em.tpr
        $GMX mdrun -deffnm ${SYSTEM_NAME}_em
        cp ${SYSTEM_NAME}.top ${SYSTEM_NAME}.itp
        sed -i -n '/\[ moleculetype \]/,$p' ${SYSTEM_NAME}.itp
        sed -i '/; Include Position restraint file/,$d' ${SYSTEM_NAME}.itp
        sed -i "s/Other_chain_A/${SYSTEM_NAME}/g" ${SYSTEM_NAME}.itp
        mkdir toppar
        cp ${SYSTEM_NAME}.itp toppar/${SYSTEM_NAME}.itp
        mkdir coord
        cp ${SYSTEM_NAME}_em.gro coord/${SYSTEM_NAME}.gro
        cp ${PETROLMD}/Quartz/files/sislab.top ${SYSTEM_NAME}.top
        sed -i "s/SiO2_2x2/${SYSTEM_NAME}/g" ${SYSTEM_NAME}.top
        $GMX editconf -f coord/${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -box 5 5 5 -noc
        $GMX solvate -cp ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top
        cp ${PETROLMD}/files/mdp-charmm36/*.mdp .
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
        $GMX mdrun -deffnm npt -update gpu
        $GMX grompp -f md_iso.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_iso.tpr
        $GMX mdrun -deffnm md_iso -update gpu

Use periodic boundary conditions
--------------------------------

One drawback of using pdb2gmx in conjunction with the speccial bonds listings is that it does not take into account periodic boundary conditions. This means that using this tooling we can not construct an ``infinite`` slab, conected to itself covalently through the periodic boundary. One way would be to add this functionality to ``pdb2gmx`` (see patch in this discussion ), but we can use the folowing trick instead.

The bonds within residue (within one unit) are added regardless of periodic boundary. Hence, if we will can move our slab molecule so that the periodic boundary goes through the molecule, instead of going through the special bonds connections. This way ``pdb2gmx`` will compute the interatomic distances correctly, thus adding the special bonds. Note that it will complain about some of the bonds being too long (these that will happen to be across the periodic boundary). It will not matter is simulations, where the periodic boundary conditions are taken into account. That is if the PBC are set correctly in the simulations.

    .. code-block:: shell

        GMX=/usr/local/gromacs/bin/gmx
        SYSTEM_NAME=slab
        FFHOME=~/git/artemzhmurov/charmm36
        PETROLMD=~/git/artemzhmurov/petrolmd
        ${PETROLMD}/build/Quartz/create_quartz
        cp ${FFHOME}/specbond.dat .
        $GMX pdb2gmx -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -ff charmm36 -water tip3p
        $GMX solvate -cp ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top
        cp ${PETROLMD}/files/mdp-charmm36/*.mdp .
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
        $GMX mdrun -deffnm npt -update gpu
        $GMX grompp -f md_anis.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_anis.tpr
        $GMX mdrun -deffnm md_anis -update gpu -nsteps 500000