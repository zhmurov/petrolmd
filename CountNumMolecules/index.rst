Simulating arbitrary mixture of molecules
=========================================

Running simulations
-------------------

To run simulations of arbitrary mixture of hydrocarbons, you are going to need:
``toppar`` folder with molecular topologies of the system with the matching coordinates in ``coord`` folder.
Both are available in the forcefield repository.
You will need to create topology and input file for the Packmol software.
This can be either done manually, or using the ``count_mols`` program.
Later takes the ``.dat`` file with the mixture composition as an input.
In this file, we list the names of the molecules and their mass concentration (in percent).
The code will count the number of molecules assuming the density of 1 kg/l based on the dimensions of the box provided in the input.
The propotions of the molecules will be computed based on their mass concentration listed in the file.
To separate molecules in the resulting simulation box, you can either edit the generated Packmol input, or by providing the sizes of the compartments for each component in the ``.dat`` file.
These data is optional and if not in the file, the software will assume that the molecules should be in the box with diagonal (0,0,0) to (Lx, Ly, Lz), with the later numbers provided in the input.
The ``.dat`` file for a 5/95 mixture of ethane with octane follows.

    .. code-block:: text

        C2H6        5.0
        C8H18       95.0

If we want to separate the gas (ethane) from the liquid (octane), we need to specify different volumes for these:

    .. code-block:: text

        C2H6        5.0      0.0   0.0   0.0   50.0  100.0  100.0
        C8H18       95.0    50.0   0.0   0.0  100.0  100.0  100.0

Note that the number of molecules will still be computed based on the input values of Lx, Ly, Lz, which allows one to adjust it.
Packmol may have hard time to fill the volume, hence we can provide smaller Lx, Ly and Lz to reduce the number of molecules.
The pressure control in simulation will handle the difference between the desired and real pressure when we start the simulations.

You may also find useful to have the configuration ``.mdp`` files for simulations.
These are provided as well.
Note that different forcefields require different setup for non-bonded forces, so use the correct ``.mdp`` files for your forcefield.

The ``count_mols`` program also requires the atomic weights database, which is used to convert the mass density to number of molecules.
There is a file with standard atomic weights provided with the script.

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