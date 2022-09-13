Creating alkanes
================

Create PDB files
----------------

.. code-block:: shell

    mkdir files
    cd files
    ../build/CreateAlkanes/create_alkanes

Copy additional files, e.g. iso-butane and iso-pentane pdbs (see ``files/PDBs`` folder in this repo):

.. code-block:: shell

    cp ../CreateAlkanes/files/PDBs/C4H10_ISO.pdb .
    cp ../CreateAlkanes/files/PDBs/C5H12_ISO.pdb .

Make topologies
---------------

.. code-block:: shell

    bash ../CreateAlkanes/create_topologies.sh

This will create ``.itp`` files for all the coordinates that we have in the folder.

.. code-block:: shell

    ../build/CountNumMolecules/count_mols ../CountNumMolecules/files/yamburg_recomb.dat ../CountNumMolecules/files/atomic_weights.dat

.. code-block:: shell

    $PACKMOL < yamburg_packmol.inp

Build and install GROMACS
=========================

1. Get the code:

    .. code-block:: shell

        git clone https://gitlab.com/gromacs/gromacs.git
        cd gromacs
        git checkout v2022.1

2. Build and install

    .. code-block:: shell
        
        mkdir build
        cmake -S. -Bbuild
        cd build
        make -j4
        sudo make install

3. You can optionally check the build before installing by running the tests:

    .. code-block:: shell
        
        make check -j4

GROMACS simulations of a box of 1000 octane molecules
=====================================================

Building the system
-------------------

1. Get the forcefield and make it available to GROMACS. Assuming that GROMACS in installed at ``/usr/local/gromacs``:

    .. code-block:: shell

        git clone https://github.com/wesbarnett/trappeua.git
        sudo cp -pr trappeua/trappeua.ff /usr/local/gromacs/share/gromacs/top/

2. Get the PDB file for octane molecule, e.g. from the output of ``reate_alkanes`` script above. You can also find coordinates online, for instance `here <https://www.angelo.edu/faculty/kboudrea/molecule_gallery/01_alkanes/00_alkanes.htm>`_.

3. Create topology file for a single molecule:
    
    .. code-block:: shell
        
        GMX=/usr/local/gromacs/bin/gmx
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

4. To create the coordinates for a box of octane molecules, we are going to use Packmol software. You will need ``gfortran``, which you can install by running ``sudo apt install gfortran``. To get and install Packmol:

    .. code-block:: shell

        git clone https://github.com/m3g/packmol.git
        cd packmol
        git checkout v20.3.5
        ./configure
        make
        PACKMOL=$(pwd)/packmol

5. To create a 10nm x 10nm x 10nm box containing 1000 octane molecules with Packmol, create `C8H18_1000.inp` file with the following:

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

6. Create topology file for GROMACS. First, run ``pdb2gmx`` to create ``.gro`` file and a stub for topology file. We are going to use the topology for a single octane molecule, but having a ``.top`` file to start with should help:

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
------------------

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


Creating topologies for isobutane and isopentane molecules
==========================================================

    .. code-block:: shell

        $GMX pdb2gmx -f C4H10_ISO.pdb -o C4H10_ISO.gro -p C4H10_ISO.top -i C4H10_ISO_posre.itp
        $GMX pdb2gmx -f C5H12_ISO.pdb -o C5H12_ISO.gro -p C5H12_ISO.top -i C5H12_ISO_posre.itp
