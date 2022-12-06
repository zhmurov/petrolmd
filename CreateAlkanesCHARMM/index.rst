Create topologies and coordinates for the alkanes in CHARMM36
-------------------------------------------------------------

Make sure to use patched CHARMM36 forcefield from here: https://gitlab.com/artemzhmurov/charmm36 . 

The following will create initial coordinates for alkanes, build the topologies and minimize coordinates using GROMACS:

    .. code-block:: shell

        mkdir toppar
        cd toppar
        ${PETROLMD}/build/CreateAlkanesCHARMM/create_alkanes_charmm
        bash ${PETROLMD}/CreateAlkanesCHARMM/create_topologies.sh
        bash ${PETROLMD}/CreateAlkanesCHARMM/minimize.sh

This will produce two folders: coord, with minimized coordinates in gro and pdb formats and toppar with itp files for the molecules.