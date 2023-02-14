Linear hydrocarbons (alkanes) 
=============================

The aim of this section is to create coordinates and molecular topologies for all alkanes that may be found in the mixture.
Likely there are some already present in the force-field and these will be used as a source of parameters and atom types.

Coordinates
-----------

Alkanes are chains of carbon atoms saturated with hydrogens.
Let us first look at the butane topology file to extract all the necessary parameters for the geometry of alkanes.
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

The first three rows here are the respective atom types, followed by the function id (5 is for urey-bradley angle with both harmonic angle and harmonic distance contribution).
Urey-Bradley angle has three parameters: the equilibrium angle (``theta0``), angle spring constant (``ktheta``), equilibrium distance between first and third atom (``ub0``) and distance spring constant (``kub``).

There is slightly more variability in angles for different atom types.
But the angles are close enough to rely on the energy minimization to fix inconsistencies if there will be any.
We also need to satisfy Urey-Bradley distances, but for simplicity we are going to leave that to energy minimization algorithm as well.

The dihedral angles for our carbon atom types are all having the equilibrium dihedral angles of 0 or 180 degrees, which indicates that the carbon backbone structure is planar.
Hence, we can build it in one plane, leaving e.g. z coordinates zero.
The hydrogens are not in plane.
The backbone hydrogens are sticking equidistantly from the plane while being apart from respective carbons in plane.
The terminal hydrogens should form a tetrahedral with the first (last) carbon and three hydrogens in the corners.

    .. code-block:: text

        [ dihedraltypes ]
        ...
        ;      i        j        k        l  func         phi0         kphi  mult
        ...
        CC32A    CC32A    CC32A    CC32A     9     0.000000     0.470742     5 ; alkane, c27r klauda et al 2004
        CC32A    CC32A    CC32A    CC32A     9     0.000000     0.395723     4 ; alkane, c27r klauda et al.2004
        CC32A    CC32A    CC32A    CC32A     9   180.000000     0.626554     3 ; alkane, c27r klauda et al 2004
        CC32A    CC32A    CC32A    CC32A     9     0.000000     0.269868     2 ; alkane, c27r klauda et al 2004

Let us start with carbon atoms, that form a sawtooth-like structure.
We are going to place the first carbon in (0,0,0).
If the x is the general direction of the chain, each next carbon is going to be   :math:`r_{CC}\times\sin(\alpha_{CCC})` further away from the starting point.
The y coordinates will be :math:`r_{CC}\times\cos(\alpha_{CCC})` for the odd atoms and zero for the even atoms, forming a sawtooth-like structure.
Here, :math:`r_{CC}` is the equilibrium length of the covalent bond between two carbons, :math:`\alpha_{CCC}` is the equilibrium angle between two such bonds.

The two hydrogens that are connected to the carbon in chain are :math:`\Delta y=r_{CH}\times\cos(\alpha_{HCH})` further away from the backbone in plane and are sticking out by :math:`\Delta z=r_{CH}\sin(\alpha_{HCH})` out of plane (in z direction).
Note that due to the geometry of the backbone, :math:`\Delta y` should be added to the y coordinate of the respective carbon for odd carbons and subtracted for the even carbons.
The :math:`\Delta z` should also be added and subtracted from the z coordinates of the respective carbon for the two connected hydrogens.

The terminal carbons are connected to three hydrogens, with four atoms forming a tetrahedron.
To simplify the geometry, we are going to assume that the axis of this tetrahedron is along the x (i.e. along the main axis of the alkane, not along the first C-C covalent bond).
In this case, the x coordinate of all hydrogens should be a height of the tetrahedron away from the carbon atom, i.e. :math:`x = r_{CH}\sqrt\left(1-\frac{4}{3}\sin^2(\alpha_{HCH}/2)\right)`.
If the vertex of tetrahedron is on the same z-plane as the carbon atom, its y coordinate differ by :math:`2r_{CH}\sin(\alpha_{HCH}/2)`, otherwise by :math:`r_{CH}\sin(\alpha_{HCH}/2)`. In the later case, difference on z axis is :math:`2r_{CH}\cos(\alpha_{HCH}/2)`.
Here, :math:`r_{CH}` is the equilibrium length of the covalent bond between carbon and hydrogen, :math:`\alpha_{HCH}` is the equilibrium angle between teo C-H bonds.

One other special case is a methane molecule, in which case the hydrogens are in the vertices of a `regular tetrahedron <https://en.wikipedia.org/wiki/Tetrahedron#:~:text=(Vertex%20figure),edges%2C%20and%20four%20vertex%20corners.>`_.

Now we can create coordinates for alkane chain of an arbitrary length, let us build topologies for these molecules.

Topology
--------

We already looked at the topology of the butane molecule.
Making topologies for the rest of alkanes should be then simple using analogy.
In fact, there are some more in the ``ethers.itp`` file: ethane, octane. decane, etc.
However, we can take it a notch further, by using the same mechanism that GROMACS uses for proteins.
Proteins are sequences of amino acid residues, linked covalently.
Since each residue is one of 20 standard amino-acids, there are :math:`20^100` variants to construct a 100-residue-long protein.
Obviously, it is not possible to keep topology for all of them.
Instead, each residue has special entries in their bond listings on how to connect the residue to the next/previous one, marked with +/- signs.
Similarly, we can construct alkanes out of `residues`: CH3 and CH2 groups.
For these, we can also add the intra-residue bonds, so that they will all be connected in the chain.
To do this, add an ``alkanes.rtp`` file to the force-field folder and add following lines into it:

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

Here we have residue ``CH2``, containing of three atoms: one carbon (``C``) and two hydrogens (``H1`` and ``H2``).
The atom types are listed in the second column of the ``atoms`` section, followed by the atom charges.
These are taken from the butane topology (see above).
The last column is the charge group: the total charge of a group should be integer number of electron charges (zero in our case).
The atom names should be unique, since they define the connectivity within and between residues.
As it is follows from the ``bonds`` section, the carbon atom is connected to both hydrogens and to carbon atom of residue before (``-C``) and after (``+C``).

It is very important that the atom names used here are matched with those we specify in the coordinate files.
Also, it is important the the residues are named and numbered properly there.
Here is an example of such file for butane in ``.gro`` format:

    .. code-block:: text

        C4H10
        14
            1CH3      C    1   0.129   0.169   0.182
            1CH3     H1    2   0.060   0.254   0.182
            1CH3     H2    3   0.108   0.108   0.273
            1CH3     H3    4   0.108   0.108   0.091
            2CH2      C    5   0.263   0.241   0.182
            2CH2     H1    6   0.266   0.307   0.272
            2CH2     H2    7   0.266   0.307   0.092
            3CH2      C    8   0.389   0.157   0.182
            3CH2     H1    9   0.387   0.091   0.273
            3CH2     H2   10   0.387   0.091   0.091
            4CH3      C   11   0.523   0.229   0.182
            4CH3     H1   12   0.593   0.145   0.182
            4CH3     H2   13   0.544   0.291   0.273
            4CH3     H3   14   0.544   0.291   0.091

There are 4 residues in total, in the following order: CH3-CH2-CH2-CH3.

Note that these topologies above don't have the angle and dihedral angle description, neither there is a type for the bonds.
These are constructed using the criteria that should be added to the beginning of the ``alkanes.rtp`` file:

    .. code-block:: text
                
        [ bondedtypes ]
        ; Col 1: Type of bond 
        ; Col 2: Type of angles 
        ; Col 3: Type of proper dihedrals 
        ; Col 4: Type of improper dihedrals 
        ; Col 5: Generate all dihedrals if 1, only heavy atoms of 0. 
        ; Col 6: Number of excluded neighbors for nonbonded interactions 
        ; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1 
        ; Col 8: Remove propers over the same bond as an improper if it is 1 
        ; bonds  angles  dihedrals  impropers  all_dihedrals  nrexcl  HH14  RemoveDih 
        1       5        9          2            1           3      1       0

As it follows from the column-by-column comments, GROMACS will generate harmonic bonds (type 1), Urey-Bradley angles (type 5), multiple proper dihedrals (type 9), etc. (more information on bond type is available in `GROMACS manual <https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html>`_).
Column 5 indicates that all dihedrals will be generated, which we will see once it is done.

Another corner-case to take care is the methane molecule, where 4 hydrogens are connected to the carbon atom.
The residue description for this molecule is:

    .. code-block:: text

        ; Methane
        [ CH4 ]
        [ atoms ]
            C      CC33A -0.3600   1
            H1     HCA3A  0.0900   1
            H2     HCA3A  0.0900   1
            H3     HCA3A  0.0900   1
            H4     HCA3A  0.0900   1
        [ bonds ]
            C    H1
            C    H2
            C    H3
            C    H4

Generate topology and run energy minimization
---------------------------------------------

Once the coordinates and ``alkanes.rst`` files are ready, we can repeat the procedure we used for the benzene molecule.
First, we create a full topology or the system using ``gmx pdb2gmx`` utility, then modify it, splitting out the molecular topology.
Then we can run energy minimization in vacuum and save energy minimized coordinates in ``.gro`` and ``.pdb`` formats.
Later is useful for third-party software, e.g. for Packmol, if we are to create a mixture of compounds.
We can also do these procedures in bulk, if we create coordinates for alkanes of different length.
The following script rely on properly set environment variables to define path to GROMACS and name of the molecule.

    .. code-block:: text

        # Run GROMACS to create top files for all the PDBs in the folder
        $GMX pdb2gmx -f ${name}.gro -o ${name}.gro -p ${name}.top -i ${name}_posre.itp -ff charmm36 -water tip3p
        # Create a copy of the topology that can be included
        cp ${name}.top ${name}.itp
        # Remove th1e header
        sed -i -n '/\[ moleculetype \]/,$p' ${name}.itp
        # Remove the footer
        sed -i '/; Include Position restraint file/,$d' ${name}.itp
        # Rename the molecule
        sed -i "s/Other/${name}/g" ${name}.itp
        # Combine topologies into one file
        # cat ${name}.itp >> alkanes.itp
        # Copy the topolgy to separate folder
        cp ${name}.itp toppar/${name}.itp
        
        # Minimize energy
        
        # Move molecule so that all coordinates are positive
        $GMX editconf -f ${name}.pdb -o ${name}.gro -d 0.1
        # Create large simulation box (~no PBC simulations)
        $GMX editconf -f ${name}.gro -o ${name}.gro -box 100 100 100 -noc
        # Configure and run GROMACS
        $GMX grompp -f em.mdp -c ${name}.gro -p ${name}.top -o ${name}_em.tpr
        $GMX mdrun -deffnm ${name}_em
        # Rename the system in the resulting file
        sed -i "s/Protein/${name}/g" ${name}_em.gro
        # Convert GRO to PDB
        $GMX editconf -f ${name}_em.gro -o ${name}_em.pdb
        # Copy the resulting coordinates
        cp ${name}_em.gro coord/${name}.gro
        cp ${name}_em.pdb coord/${name}.pdb

Create topologies and coordinates for the alkanes in CHARMM36
-------------------------------------------------------------

Make sure to use patched CHARMM36 forcefield from here: https://gitlab.com/artemzhmurov/charmm36 , or you patched the forcefield as it is described above (also see version of the force-field provided here).
The script will use the code that create coordinate files for alkanes starting from CH4 up to C100H202 (you can change it in the code, see ``main.cpp`` in this folder).
There are also pre-generated ``.pdb`` files for isobutane and isopentane, that require additional patch to the forcefield:

    .. code-block:: text

        ; Methyl branch connected to CH
        [ MB ]
        [ atoms ]
            C1      CC33A -0.2700   1
            H11     HCA3A  0.0900   1
            H12     HCA3A  0.0900   1
            H13     HCA3A  0.0900   1
            C       CC32A -0.0900   1
            H1      HCA2A  0.0900   1
        [ bonds ]
        C1    H11
        C1    H12
        C1    H13
        C     C1
        C     H1
        -C     C
        C     C+

The following script will create initial coordinates for alkanes, build the topologies and minimize coordinates using GROMACS.
Make sure that you change the path to GROMACS and to this repo in the script before running.

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
