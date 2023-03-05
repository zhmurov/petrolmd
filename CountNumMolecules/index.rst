Simulating arbitrary mixture of molecules
=========================================

Running simulations
-------------------

To run simulations of arbitrary mixture of hydrocarbons, you are going to need:
A folder with molecular topologies of the system with the matching coordinates (available in the forcefield repository).
You will need to create topology and input file for the Packmol software.
This can be either done manually, or using the ``count_mols`` program.
Later takes the ``.dat`` file with the mixture composition as an input.
In this file, we list the names of the molecules and their mass concentration (in percent).
The code will count the number of molecules assuming the density of 1 kg/l based on the dimensions of the box provided in the input.
The proportions of the molecules will be computed based on their mass concentration listed in the file.
To separate molecules in the resulting simulation box, you can either edit the generated Packmol input, or by providing the sizes of the compartments for each component in the ``.dat`` file.
These data is optional and if not in the file, the software will assume that the molecules should be in the box with diagonal (0,0,0) to (Lx, Ly, Lz), with the later numbers provided in the input.
The ``.dat`` file for a 5/95 mixture of ethane with octane follows.

    .. code-block:: text

        C2H6        5.0
        C8H18       95.0

To get topology and Packmol input file, run the ``count_mols`` script.

    .. code-block:: shell

        ${PETROLMD}/build/CountNumMolecules/count_mols ${PETROLMD}/CountNumMolecules/files/atomic_weights.dat ${PETROLMD}/CountNumMolecules/files/${SYSTEM_NAME}.dat ${SYSTEM_NAME} ${Lx} ${Ly} ${Lz}

These will generate the following two files:

    .. code-block:: text

        tolerance 2.0
        filetype pdb
        output ethane-octane.pdb

        structure charmm36.ff/hydrocarbons/C2H6.pdb
        number 8011
        inside box 0 0 0 200 200 200
        end structure
        structure charmm36.ff/hydrocarbons/C8H18.pdb
        number 40066
        inside box 0 0 0 200 200 200
        end structure

    .. code-block:: text

        ; Include forcefield parameters
        #include "charmm36.ff/forcefield.itp"
        #include "charmm36.ff/tip3p.itp"
        #include "charmm36.ff/hydrocarbons/C2H6.itp"
        #include "charmm36.ff/hydrocarbons/C8H18.itp"

        [ system ]
        ; Name
        ethane-octane

        [ molecules ]
        ; Compound     #mols
        C2H6      8011
        C8H18      40066



If we want to separate the gas (ethane) from the liquid (octane), we need to specify different volumes for these:

    .. code-block:: text

        C2H6        5.0      0.0   0.0   0.0   5.0  10.0  10.0
        C8H18       95.0     5.0   0.0   0.0  10.0  10.0  10.0

Note that the number of molecules will still be computed based on the input values of Lx, Ly, Lz, which allows one to adjust it.
Packmol may have hard time to fill the volume, hence we can provide smaller Lx, Ly and Lz to reduce the number of molecules.
The pressure control in simulation will handle the difference between the desired and real pressure when we start the simulations.
If using the ``.dat`` file above, the resulting ``.top`` file will be the same, Packmol script will have new boundaries for the respective mixture components.

    .. code-block:: text

        tolerance 2.0
        filetype pdb
        output ethane-octane.pdb

        structure charmm36.ff/hydrocarbons/C2H6.pdb
        number 8011
        inside box 0 0 0 50 100 100
        end structure
        structure charmm36.ff/hydrocarbons/C8H18.pdb
        number 40066
        inside box 50 0 0 100 100 100
        end structure


You may also find useful to have the configuration ``.mdp`` files for simulations.
These are provided as well.
Note that different forcefields require different setup for non-bonded forces, so use the correct ``.mdp`` files for your forcefield.

The ``count_mols`` program also requires the atomic weights database, which is used to convert the mass density to number of molecules.
There is a file with standard atomic weights provided with the script.

    .. code-block:: shell


        GMX=/usr/local/gromacs/bin/gmx
        PACKMOL=~/git/external/packmol/packmol
        PETROLMD=~/git/artemzhmurov/petrolmd/
        CHARMM36_HOME=~/git/artemzhmurov/charmm36
        Lx=20
        Ly=20
        Lz=20
        SYSTEM_NAME=yamburg_recomb

        cp -r ${CHARMM36_HOME}/toppar/ .
        cp -r ${CHARMM36_HOME}/coord/ .
        cp ${PETROLMD}/files/mdp-charmm36/*.mdp .
        ${PETROLMD}/build/CountNumMolecules/count_mols ${PETROLMD}/CountNumMolecules/files/atomic_weights.dat ${PETROLMD}/CountNumMolecules/files/${SYSTEM_NAME}.dat ${SYSTEM_NAME} ${Lx} ${Ly} ${Lz}
        $PACKMOL < ${SYSTEM_NAME}_packmol.inp
        $GMX editconf -f ${SYSTEM_NAME}.pdb -o ${SYSTEM_NAME}_box.gro -box ${Lx} ${Ly} ${Lz} -noc
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}_box.gro -p ${SYSTEM_NAME}.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p ${SYSTEM_NAME}.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p ${SYSTEM_NAME}.top -o npt.tpr
        $GMX mdrun -deffnm npt -update gpu
        $GMX grompp -f md_iso.mdp -c npt.gro -p ${SYSTEM_NAME}.top -o md_iso.tpr
        $GMX mdrun -deffnm md_iso -update gpu