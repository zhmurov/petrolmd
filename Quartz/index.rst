Creating the coordinates and topology for a quartz slab of a given size
-----------------------------------------------------------------------

    .. code-block:: shell

        cp ~/git/artemzhmurov/charmm36/specbond.dat .
        cp ~/git/artemzhmurov/petrolmd/files/em_vac.mdp em.mdp
        ~/git/artemzhmurov/petrolmd/build/Quartz/create_quartz
        $GMX pdb2gmx -f slab.pdb -o slab.gro -p slab.top -ff charmm36 -water tip3p
        $GMX editconf -f slab.gro -o slab.gro -d 0.1
        $GMX editconf -f slab.gro -o slab.gro -box 100 100 100 -noc
        $GMX grompp -f em.mdp -c slab.gro -p slab.top -o slab_em.tpr
        $GMX mdrun -deffnm slab_em
        cp slab.top slab.itp
        sed -i -n '/\[ moleculetype \]/,$p' slab.itp
        sed -i '/; Include Position restraint file/,$d' slab.itp
        sed -i "s/Other_chain_A/SiO2_2x2/g" slab.itp
        mkdir toppar
        cp slab.itp toppar/SiO2_2x2.itp
        cp slab_em.gro sislab.gro
        $GMX editconf -f sislab.gro -o sislab.gro -box 5 5 5 -noc
        $GMX solvate -cp sislab.gro -o sislab_solv.gro -p sislab.top
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

    .. code-block:: shell

        ~/git/artemzhmurov/petrolmd/build/Quartz/create_quartz
        cp ~/git/artemzhmurov/charmm36/specbond.dat .
        $GMX pdb2gmx -f slab.gro -o slab.gro -p slab.top -ff charmm36 -water tip3p
        $GMX solvate -cp slab.gro -o slab_solv.gro -p slab.top
        cp ${PETROLMD}/files/mdp-charmm36/*.mdp .
        $GMX grompp -f em.mdp -c slab_solv.gro -p slab.top -o em.tpr
        $GMX mdrun -deffnm em
        $GMX grompp -f nvt.mdp -c em.gro -p slab.top -o nvt.tpr
        $GMX mdrun -deffnm nvt
        $GMX grompp -f npt.mdp -c nvt.gro -p slab.top -o npt.tpr
        $GMX mdrun -deffnm npt -update gpu
        $GMX grompp -f md_anis.mdp -c npt.gro -p slab.top -o md_anis.tpr
        $GMX mdrun -deffnm md_anis -update gpu -nsteps 500000