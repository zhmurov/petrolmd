GROMACS simulations of a box of 1000 octane molecules
=====================================================

System coordinates
------------------

To create the system, we will use Packmol software.
Packmol takes coordinates of a single molecule and fills the specified volume with translated and rotated copies of the molecule.
We assume that there is a ``coord`` folder with coordinates of single molecules that can be used to create a mixture.
To create a 10nm x 10nm x 10nm box containing 1000 octane molecules with Packmol, create ``C8H18_1000.inp`` file with the following:

    .. code-block:: text

        tolerance 2.0
        filetype pdb
        output C8H18_1000.pdb

        structure coord/C8H18.pdb
        number 1000 
        inside box 0. 0. 0. 100. 100. 100. 
        end structure

The first three lines specify the output parameters.
Distance tolerance (in angstroms) is the minimal distance between molecules.
Packmol will try to put molecules in the specified dimensions so that the minimal distance between two molecules is larger than tolerance.
Next two lines specify the output format and file name.
For each component, we specify the input file, number of molecules to add and the dimensions of the box to put molecules in (in angstroms).
To create the PDB file with Packmol, feed the created file to the software:

    .. code-block:: shell
    
        $PACKMOL < C8H18_1000.inp

This should create a ``C8H18_1000.pdb`` file. Feel free to load it into VMD or other visualization software to have a look.

System topology
---------------

Now we have coordinates for the system, we need to create topology file.
Since we already created the topology files for all the molecules in the mixture (octane in this case), the system topology is quite simple.
All we need to do is to include pre-generated topologies for molecules in the mixture and list the number of molecules.
Note that the number of molecules and their quantities should be the same as used when the coordinates were generated. The resulting top file is:

    .. code-block:: text

        ; Include forcefield parameters
        #include "charmm36.ff/forcefield.itp"
        #include "charmm36.ff/tip3p.itp"
        #include "toppar/C8H18.itp"

        [ system ]
        ; Name
        1000 octane molecules

        [ molecules ]
        ; Compound        #mols
        C8H18             1000

Adding water solvent
--------------------

There are two ways to add water molecules to the system: using GROMACS ``solvate`` utility or using Packmol with water as another component of the mixture.
To add water using GROMACS, create ``.gro`` file with solvation box specified and run ``solvate`` utility:

    .. code-block:: shell
    
        $GMX editconf -f C8H18_1000.pdb -o C8H18_1000_box.gro -box 10 10 10
        $GMX solvate -cp C8H18_1000_box.gro -cs tip4p.gro -o C8H18_1000_solv.gro -p C8H18_1000.top

    This will specify the box size to 10nm x 10nm x 10nm and add as many water molecules as it will be able to without introducing steric clashes between molecules and without exceeding normal density.
    You can also specify maximum number of water molecules by adding ``--maxsol`` parameter to ``gmx solvate``.
    See ``gmx solvate`` `section in GROMACS manual <https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html>`_ for more details on how to use this utility.
    Note that this will overwrite the ``.top`` file, adding the solvent (water) molecules into the system description.
    You can edit the name of the system and/or rename the final topology file if you wish.

Another way of adding water molecules is to deal with them the same way as with other components of the mixture with Packmol.
To do so, you will need a coordinate file for single water molecule, which can be created similarly to the rest of the components.
To mix 1000 octane molecules with 5000 water molecules, one can use the following Packmol input file:

    .. code-block:: text

        tolerance 2.0
        filetype pdb
        output C8H18_1000_solv.pdb

        structure coord/C8H18.pdb
        number 1000 
        inside box 0. 0. 0. 100. 100. 100. 
        end structure

        structure coord/tip3p.pdb
        number 5000 
        inside box 0. 0. 0. 100. 100. 100. 
        end structure

The corresponding topology file will be:

    .. code-block:: text

        ; Include forcefield parameters
        #include "charmm36.ff/forcefield.itp"
        #include "charmm36.ff/tip3p.itp"
        #include "toppar/C8H18.itp"

        [ system ]
        ; Name
        1000 octane molecules in water

        [ molecules ]
        ; Compound        #mols
        C8H18             1000
        SOL               5000

The drawback of this approach is that user has to pre-compute the number of water molecules instead of relying on GROMACS to add water up to desired density.
The first advantage is that we eliminate one step in the procedure of system generation.
But more importantly, we can now easily create system where the two components are separated at the start of the simulations by providing different compartments to Packmol:

    .. code-block:: text

        tolerance 2.0
        filetype pdb
        output C8H18_1000_solv.pdb

        structure coord/C8H18.pdb
        number 1000 
        inside box 0. 0. 0. 50. 50. 50. 
        end structure

        structure coord/tip3p.pdb
        number 5000 
        inside box 50. 50. 50. 100. 100. 100. 
        end structure

This way, the simulation box will be split in half with the left side filled with octane and right side filled with water.
Note that Packmol provides an extensive number of options to specify the geometry of the system.
See `Packmol users guide <https://m3g.github.io/packmol/userguide.shtml>`_ for details.

Preparing the system
--------------------


Energy minimization

    .. code-block:: shell

        $GMX grompp -f em.mdp -c C8H18_1000_solv.gro -p C8H18_1000.top -o em.tpr
        $GMX mdrun -deffnm em

Equilibration

    NVT:

    .. code-block:: shell

        $GMX grompp -f nvt.mdp -c em.gro -p C8H18_1000.top -o nvt.tpr
        $GMX mdrun -deffnm nvt

    NPT:

    .. code-block:: shell

        $GMX grompp -f npt.mdp -c nvt.gro -p C8H18_1000.top -o npt.tpr
        $GMX mdrun -deffnm npt


Production simulations
--------------------

    .. code-block:: shell

        $GMX grompp -f md.mdp -c npt.gro -p C8H18_1000.top -o md.tpr
        $GMX mdrun -deffnm md

