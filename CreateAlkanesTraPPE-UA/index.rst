GROMACS simulations of a box of 1000 octane molecules
-----------------------------------------------------

Building the system
^^^^^^^^^^^^^^^^^^^

1. Get the PDB file for octane molecule, e.g. from the output of ``create_alkanes`` script above. You can also find coordinates online, for instance `here <https://www.angelo.edu/faculty/kboudrea/molecule_gallery/01_alkanes/00_alkanes.htm>`_.

2. Create topology file for a single molecule:
    
    .. code-block:: shell
        
        $GMX pdb2gmx -f C8H18.pdb -o C8H18.gro -p C8H18.top -i C8H18_posre.itp

    Select TraPPE-UA forcefield for hydrocarbons ("Transferable Potentials for Phase Equilibria - United Atom (TraPPE-UA) with HH-Alkane modifications") and TIP4P ("TIP4P2005  TIP4P/2005") for water.

    Create a copy of the generated ``.top`` file, named ``C8H18.itp``. This will be lated used as a component of the bigger system. In this file, remove the ``[ pairs ]`` section, remove the system description as well as all references to water model (i.e. remove everything starting from ``; Include Position restraint file``), replace the name ``Other`` with ``C8H18``:

    .. code-block:: shell

        cp C8H18.top C8H18.itp
        sed -i -n '/\[ moleculetype \]/,$p' ${name}.itp
        # Remove pairs section (needed by TraPPE forcefield)
        sed -i -n '1,/pairs/p;/angles/,$p' ${name}.itp
        sed -i '\[ pairs \]/d' C8H18.itp
        # Remove the footer
        sed -i '/; Include Position restraint file/,$d' ${name}.itp
        # Rename the molecule
        sed -i "s/Other/${name}/g" ${name}.itp

3. To create a 10nm x 10nm x 10nm box containing 1000 octane molecules with Packmol, create `C8H18_1000.inp` file with the following:

    .. code-block:: text

        tolerance 2.0
        filetype pdb
        output C8H18_1000.pdb

        structure C8H18.pdb
        number 1000 
        inside box 0. 0. 0. 100. 100. 100. 
        end structure

    And feed it to packmol executable:

    .. code-block:: shell
    
        $PACKMOL < C8H18_1000.inp

    This should create a ``C8H18_1000.pdb`` file. Feel free to load it into VMD or other visualization software to have a look.

4. Create topology file for GROMACS. First, run ``pdb2gmx`` to create ``.gro`` file and a stub for topology file. We are going to use the topology for a single octane molecule, but having a ``.top`` file to start with should help:

    .. code-block:: shell
        
        $GMX pdb2gmx -f C8H18_1000.pdb -o C8H18_1000.gro -p C8H18_1000.top -i C8H18_1000_posre.itp

    Use TraPPE-UA and TIP4P forcefields.

    We are going to use the topology of a single octane molecule that we created earlier. Hence, we don't need the description of the molecules in the topology file. 
    
    
    So, remove the lines starting from ``[ moleculetype ]`` all the way to ``[ system ]``. You can do this manually, or by executing the following command:

    .. code-block:: shell

        sed -i '/^\[ moleculetype \]/,/\[ system \]/{/\[ system \]/b;d;}' C8H18_1000.top

    Next, include the ``.itp`` file for the single octane molecule by adding:

    .. code-block:: text

        #include "C8H18.itp"

    Name the molecule appropriately and modify the description of the system to include 1000 molecules of ``C8H18``. The resulting top file should look something like this:

    .. code-block:: text

        ; Include forcefield parameters
        #include "trappeua.ff/forcefield.itp"
        #include "trappeua.ff/tip4p2005.itp"
        #include "C8H18.itp"

        [ system ]
        ; Name
        1000 octane molecules

        [ molecules ]
        ; Compound        #mols
        C8H18             1000

System preparation
^^^^^^^^^^^^^^^^^^

1. Solvate the system in water:

    .. code-block:: shell
    
        $GMX editconf -f C8H18_1000.gro -o C8H18_1000_box.gro -box 10 10 10
        $GMX solvate -cp C8H18_1000_box.gro -cs tip4p.gro -o C8H18_1000_solv.gro -p C8H18_1000.top

    Note that this will overwrite the ``.top`` file, adding the solvent (water) molecules into the system description. You can edit the name of the system and/or rename the final topology file if you wish.


2. Energy minimization

    .. code-block:: shell

        $GMX grompp -f em.mdp -c C8H18_1000_solv.gro -p C8H18_1000.top -o em.tpr
        $GMX mdrun -deffnm em

3. Equilibration

    NVT:

    .. code-block:: shell

        $GMX grompp -f nvt.mdp -c em.gro -p C8H18_1000.top -o nvt.tpr
        $GMX mdrun -deffnm nvt

    NPT:

    .. code-block:: shell

        $GMX grompp -f npt.mdp -c nvt.gro -p C8H18_1000.top -o npt.tpr
        $GMX mdrun -deffnm npt

4. Production run:

    .. code-block:: shell

        $GMX grompp -f md.mdp -c npt.gro -p C8H18_1000.top -o md.tpr
        $GMX mdrun -deffnm md


Creating alkanes
----------------

Building the helper codes
^^^^^^^^^^^^^^^^^^^^^^^^^

    .. code-block:: shell

        git clone git@gitlab.com:artemzhmurov/petrolmd.git
        cd petrolmd
        cmake -S. -Bbuild
        cmake --build build
        PETROLMD=${pwd}

Create PDB files
^^^^^^^^^^^^^^^^

    .. code-block:: shell

        mkdir toppar
        cd toppar
        ${PETROLMD}/build/CreateAlkanesTraPPE-UA/create_alkanes

Copy additional files, e.g. iso-butane and iso-pentane pdbs (see ``files/PDBs`` folder in this repo):

    .. code-block:: shell

        cp ${PETROLMD}/CreateAlkanesTraPPE-UA/files/PDBs/C4H10_ISO.pdb .
        cp ${PETROLMD}/CreateAlkanesTraPPE-UA/files/PDBs/C5H12_ISO.pdb .

Make topologies
^^^^^^^^^^^^^^^

    .. code-block:: shell

        bash ${PETROLMD}/CreateAlkanesTraPPE-UA/create_topologies.sh

This will create ``.itp`` files for all the coordinates that we have in the folder. You will also need coordinates for the water molecule:

    .. code-block:: shell

        cp ${PETROLMD}/files/tip4p.gro .
        cd ..

It is convenient to save the system name and box sizes into variables, so scripts below can be copy-pasted:

    .. code-block:: shell

        Lx=<Lx(nm)>
        Ly=<Ly(nm)>
        Lz=<Lz(nm)>
        SYSTEM_NAME=<system_name>

    .. code-block:: shell

        ${PETROLMD}/build/CountNumMolecules/count_mols ${PETROLMD}/CountNumMolecules/files/atomic_weights.dat ${PETROLMD}/CountNumMolecules/files/<composition_data>.dat ${SYSTEM_NAME} ${Lx} ${Ly} ${Lz}

This will produce two files: topology for GROMACS and input file for packmol. To create coordinates file, use:

    .. code-block:: shell

        $PACKMOL < ${SYSTEM_NAME}_packmol.inp

You should be good to go for GROMACS simulation. You can use provided ``.mdp`` files for energy minimization. equilibration and production runs:

    .. code-block:: shell

        cp ${PETROLMD}/files/*.mdp .

    .. code-block:: shell

        $GMX editconf -f ${SYSTEM_NAME}.pdb -o ${SYSTEM_NAME}_box.gro -box ${Lx} ${Ly} ${Lz}
        $GMX solvate -cp ${SYSTEM_NAME}_box.gro -cs toppar/tip4p.gro -o ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
        $GMX mdrun -deffnm npt
        $GMX grompp -f md.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md.tpr
        $GMX mdrun -deffnm md

Example script:

    .. code-block:: shell

        ${PETROLMD}/build/CountNumMolecules/count_mols ${PETROLMD}/CountNumMolecules/files/atomic_weights.dat ${PETROLMD}/CountNumMolecules/files/methane-octane.dat methane-octane 10.0 10.0 10.0
        $PACKMOL < methane-octane_packmol.inp
        $GMX editconf -f methane-octane.pdb -o methane-octane_box.gro -box 10 10 10
        $GMX solvate -cp methane-octane_box.gro -cs toppar/tip4p.gro -o methane-octane_solv.gro -p methane-octane.top
        $GMX grompp -f em.mdp -c methane-octane_solv.gro -p methane-octane.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p methane-octane.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p methane-octane.top -o npt.tpr
        $GMX mdrun -deffnm npt
        $GMX grompp -f md.mdp -c npt.gro -p methane-octane.top -o md.tpr
        $GMX mdrun -deffnm md



Creating topologies for isobutane and isopentane molecules
----------------------------------------------------------

    .. code-block:: shell

        $GMX pdb2gmx -f C4H10_ISO.pdb -o C4H10_ISO.gro -p C4H10_ISO.top -i C4H10_ISO_posre.itp
        $GMX pdb2gmx -f C5H12_ISO.pdb -o C5H12_ISO.gro -p C5H12_ISO.top -i C5H12_ISO_posre.itp


Building the system with separate compartments for water, liquid and gas phase hydrocarbons
-------------------------------------------------------------------------------------------

    .. code-block:: shell

        mkdir toppar
        cd toppar
        ${PETROLMD}/build/CreateAlkanesTraPPE-UA/create_alkanes
        cp ${PETROLMD}/CreateAlkanesTraPPE-UA/files/PDBs/C4H10_ISO.pdb .
        cp ${PETROLMD}/CreateAlkanesTraPPE-UA/files/PDBs/C5H12_ISO.pdb .
        bash ${PETROLMD}/CreateAlkanesTraPPE-UA/create_topologies.sh
        cp ${PETROLMD}/files/tip4p.gro .
        cp ${PETROLMD}/files/SOL.itp .
        cd ..
        Lx=20
        Ly=20
        Lz=20
        SYSTEM_NAME=yamburg_recomb
        $GMX solvate -cs toppar/tip4p.gro -box 10.0 20.0 10.0 -maxsol 60000
        ${PETROLMD}/build/CountNumMolecules/count_mols ${PETROLMD}/CountNumMolecules/files/atomic_weights.dat ${PETROLMD}/CountNumMolecules/files/${SYSTEM_NAME}.dat ${SYSTEM_NAME} 10 20 10

At this stage, we need to manually edit the configuration file for packmol to separate the gas phase from the liquid phase of hydrocarbons. To do so, we will edit the boxes in which the packmol will be placing the molecules in the packmol config file, generated by the ``count_mols`` script. The gas phase is normally up to the C5H12. This can be placed in (0, 0, 100; 100, 200, 200), the rest can be left in the (100, 0, 0; 200, 200, 200) box. This can be easily done by using replace-all option of your text editor. Note that packmol uses angstroms, not nanometers.

    .. code-block:: shell

        $PACKMOL < ${SYSTEM_NAME}_packmol.inp
        cp ${PETROLMD}/files/mdp-trappeua/*.mdp .
        $GMX editconf -f ${SYSTEM_NAME}.pdb -o ${SYSTEM_NAME}_box.gro -box ${Lx} ${Ly} ${Lz} -noc
        $GMX solvate -cp ${SYSTEM_NAME}_box.gro -cs toppar/tip4p.gro -o ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top -maxsol 60000
        
Now we need to transfer coordinates of the water molecules from the pre-generated water box into our system. The former are in the ``out.gro`` file, the later - in the ``${SYSTEM_NAME}_solv.gro`` file. Since we imposed the limit on the number of solvent molecules in both cases and assuming that this number was reached, it should be the same in both cases.

    .. code-block:: shell

        top2psf ${SYSTEM_NAME}.top toppar/ ${SYSTEM_NAME}.psf
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}_solv.gro -p ${SYSTEM_NAME}.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
        $GMX mdrun -deffnm npt