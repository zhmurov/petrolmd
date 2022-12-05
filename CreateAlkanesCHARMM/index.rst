Create topologies and coordinates for the alkanes in CHARMM36
-------------------------------------------------------------

Make sure to use patched CHARMM36 forcefield from here: https://gitlab.com/artemzhmurov/charmm36 . 

The following will create initial coordinates for alkanes, build the topologies and minimize coordinates using GROMACS:

    .. code-block:: shell

        mkdir toppar
        cd toppar
        ${PETROLMD}/build/CreateAlkanesCHARMM/create_alkanes_charmm
        ${PETROLMD}/CreateAlkanesCHARMM/create_topologies.sh
        ${PETROLMD}/CreateAlkanesCHARMM/minimize.sh

