Create topology and coordinates for benzene molecule
====================================================

For the molecular system that contains a mixture of small molecules it is convenient to have topology and coordinate files for each component.
This way the system coordinates can be created by a third-party software (e.g. Packmol or VMD) and the topology file will just list the molecules and their quantities, not the interatomic connectivity in the molecules.
We will use an example of benzene to explain the procedure of creating a pair of coordinates and topology files for a single molecule that can be later used when the molecule is added to the mixture of many molecules.

The coordinates
---------------

We will write basic code to synthesize the coordinates according to the atomic parameters.
Another approach would be to find the coordinates on the internet, as a stand-alone file or part of another molecular system.
The later for benzene may be, e.g. a part of phenylalanine side-chain form any protein that contains this amino-acid.
Note that the the carbon, that connects the side-chain to the backbone of the protein will have to be substituted with a hydrogen.
Also, if the system is synthesized with the correct force-field parameters, there is no need to minimize the energy (but we will do it to show the complete procedure).

Let us use CHARMM 36 force-field for the reference. If we find benzene molecule there (`BENZ` entry in cgenff.rtp file), we see that all the carbons have type `CG2R61`, hydrogens are of type `HGR61`:

    .. code-block:: text

        [ BENZ ]
        ; C6H6, benzene, adm jr.
        [ atoms ]
            CG   CG2R61 -0.1150   1
            HG    HGR61  0.1150   1
            CD1   CG2R61 -0.1150   1
            HD1    HGR61  0.1150   1
            CD2   CG2R61 -0.1150   1
            HD2    HGR61  0.1150   1
            CE1   CG2R61 -0.1150   1
            HE1    HGR61  0.1150   1
            CE2   CG2R61 -0.1150   1
            HE2    HGR61  0.1150   1
            CZ   CG2R61 -0.1150   1
            HZ    HGR61  0.1150   1
        [ bonds ]
            CD1    CG
            CD2    CG
            CE1   CD1
            CE2   CD2
            CZ   CE1
            CZ   CE2
            CG    HG
            CD1   HD1
            CD2   HD2
            CE1   HE1
            CE2   HE2
            CZ    HZ


These atom types are now can be used to extract equilibrium parameters for bonds and angles form `ffbonded.itp` file in the force-field. Note that we are licky that our molecule is planar and hence we do not have to think on satisfying the dihedral angles. The parameters that we are looking for are in the following lines of the `ffbonded.itp` file:

    .. code-block:: text
        
        ...
        CG2R61   CG2R61     1   0.13750000    255224.00 ; PROT benzene, JES 8/25/89
        ...
        CG2R61    HGR61     1   0.10800000    284512.00 ; PROT phe,tyr JES 8/25/89
        ...
        CG2R61   CG2R61   CG2R61     5   120.000000   334.720000   0.24162000     29288.00 ; PROT JES 8/25/89
        ...
        CG2R61   CG2R61    HGR61     5   120.000000   251.040000   0.21525000     18409.60 ; PROT JES 8/25/89 benzene
        ...

As one can see from the comments, these atom types are found in phenylalanine and tyrosine amino-acids as well as in benzene molecule.
The equilibrium bond distances are 0.1375 nm for bond between carbons and 0.108 for carbon-hydrogen bond. The equilibrium angle is 120 degrees for both C-C-C and C-C-H angles.

So we need to create a regular hexagon with carbons in its corners and hydrogens sticking out.
Each side of the hexagon is a=0.1375 nm and hydrogens are farther b=0.108 nm away from the corners.

If placed in the x-y plane so that the center of hexagon is (0,0,0) and the longest diagonal along x-axis, the coordinates of the vertices are (clockwise):

    .. math::

        (a,0,0), (a/2,-cos(30)*a,0), (-a/2,-cos(30)*a,0), (-a,0,0), (-a/2,cos(30)*a,0), (a/2,cos(30)*a,0)

The coordinates of the respected hydrogens are in the same layout, but (a+b) away from the center, i.e.:

    .. math::

        (a+b,0,0), ((a+b)/2,-cos(30)*(a+b),0), (-(a+b)/2,-cos(30)*(a+b),0), (-(a+b),0,0), (-(a+b)/2,cos(30)*(a+b),0), ((a+b)/2,cos(30)*(a+b),0)
        
We also need to synchronize the coordinates with the topology: molecular simulation software needs to know which atom is connected to which.
This is described in the `[ bonds ]` section of the topology: atom `CD1` is connected to atom `CG`, `CD2` to `CG` as well, and so on.
Hence, if we want to go around the carbon ring, starting from the `CG` atom, the atom order will be `CG`, `CD1`, `CE1`, `CZ`, `CE2`, `CD2`.
This notation is inherited from proteins, where atoms are named according to their separation from the C-alpha atom: C-beta (CB), C-gamma (CG), C-delta (CD) and so on.
If the chain separates, atom names receive a number (e.g. `CD1` and `CD2`).
This does not make much sense in case of benzene molecule, but this does not affect us in the future so we can leave it as it is.

Once we decide which coordinates correspond to which atom name, we can create the coordinates file in any format that is suitable for the software.
For instance, we can use `GROMACS .gro format <https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#gro>`_.


Molecular topology
------------------

Once the coordinates are ready, the topology file can be created.
If the atoms and residue are named according to the force-field, the basic topology can be created using pdb2gmx utility from GROMACS:

    .. code-block:: bash

        GMX=/usr/local/gromacs/bin/gmx
        SYSTEM_NAME=C6H6
        $GMX pdb2gmx -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -ff charmm36 -water tip3p

We set and use some variable to make it easier to copy-paste these commands later.
Here, we take the .gro file we created and pass it to the `pdb2gmx`.
We use charmm36 force-field and tip3p water model (if we don't specify these options, they will be asked interactively).
As an output we have corrected .gro file (if the correction is needed, e.g. hydrogen were added) and topology file in .top format.
We opted to overwrite the .gro file here, since not much is updated there --- only the periodic boundary description is added.
Don't worry if you wanted to compare the files --- GROMACS makes back-ups it case it overwrites (see the log file for the backup name).

Before we modify the topology to make it more useful in the future, we can finalize our coordinates by running energy minimization simulations.
Strictly speaking, this step is not required because we created the coordinates using force-field parameters.
Hence the structure should be already in its energy minima.
But it is useful to do the energy minimization to fix any imperfections in the structure and to make sure that the coordinates are useable.

The following script performs the energy minimization in vacuo for the molecular system.
We first run the `editconf` utility twice to set the periodic boundary correctly.
The first run moves the molecule so that all the coordinates are positive.
The `-d 0.1` option means, that the box will be no closer than 0.1 nm from the molecule.
By default, the `editconf` shifts the coordinates to the center of the box, which is in the positive quadrant.
The second run of `editconf` changes the periodic box definition to a cube with 100 nm side.
This is an arbitrary large number to make sure that PBC don't affect the system --- modern versions of GROMACS can only run with periodic boundary.
Note that we use `-noc` option here to leave the molecule where it is instead of moving it to the center of the box.
After coordinates are prepared, we configure GROMACS with `grompp` utility.
This takes the .mdp file that describes the simulation protocol, topology and coordinate files.
`grompp` creates a portable `.tpr` file, which contains all the data for the molecular dynamics simulation run.
This file is used with `gmx mdrun`, which is the main simulation engine.
The `-deffnm` option means default name, i.e. the name of the files for input/output that will be different only by extension.
Hence, the resulting (energy minimized) structure will be saved as `${SYSTEM_NAME}_em.gro`, which are the coordinates that we are going to save for future use.

    .. code-block:: bash

        $GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -d 0.1
        $GMX editconf -f ${SYSTEM_NAME}.gro -o ${SYSTEM_NAME}.gro -box 100 100 100 -noc
        cp ${PETROLMD}/files/em_vac.mdp em.mdp
        $GMX grompp -f em.mdp -c ${SYSTEM_NAME}.gro -p ${SYSTEM_NAME}.top -o ${SYSTEM_NAME}_em.tpr
        $GMX mdrun -deffnm c
        mkdir coord
        cp ${SYSTEM_NAME}_em.gro coord/${SYSTEM_NAME}.gro

Current `.top` file does two things at the same time: it describes the composition of the entire molecular system and topology of a single benzene molecule.
This is not surprising because our system is singular benzene molecule.
However, this is not convenient if we want to re-use the molecular topology description part of the file when benzene molecule(s) is(are) part of another molecular system.
What we can do in this case is we can create a separate `.itp` file and copy the topology description there.
The part we need to keep starts with `[ moleculetype ]` and ends after `[ dihedrals ]` section of the topology file.
By default, the molecule is named `Other`, which can also be changed for convenience (and should be change to avoid duplicating names).
The described procedures can be performed with the following script:

    .. code-block:: bash

        cp ${SYSTEM_NAME}.top ${SYSTEM_NAME}.itp
        sed -i -n '/\[ moleculetype \]/,$p' ${SYSTEM_NAME}.itp
        sed -i '/; Include Position restraint file/,$d' ${SYSTEM_NAME}.itp
        sed -i "s/Other/${SYSTEM_NAME}/g" ${SYSTEM_NAME}.itp
        mkdir toppar
        cp ${SYSTEM_NAME}.itp toppar/${SYSTEM_NAME}.itp

Now, the system topology can just include the molecular topology file.
For the system of a singular benzene, the topology will be:

    .. code-block:: text

        ; Include forcefield parameters
        #include "charmm36.ff/forcefield.itp"
        #include "toppar/C6H6.itp"

        [ system ]
        ; Name
        C6H6

        [ molecules ]
        ; Compound     #mols
        C6H6      1

Creating topologies and coordinates for all the molecules of the mixture allow to re-use them in the future, thus simplifying the process of simulating the arbitrary mixture significantly.
The coordinates file can be constructed using pre-generated coordinate files for the molecules (e.g. using Packmol software) and the topologies of molecules should be included into the system topology file.
These two files should be kept in-sync though, with the number of molecules and their order in coordinates file strictly correspondent to the number of molecules and their order in the system topology.

The .mdp file
-------------

Above, we ran something called energy minimization in vacuum.
This is not something that is done by default, we explicitly asked molecular simulation software to do it.
This is specified in the configuration or molecular dynamics parameters (`.mdp`) file (see full description `here <https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html>`_).
We used the following parameters in our run:

    .. code-block:: text

        title       = enrgy minimisation

        ; Parameters describing what to do, when to stop and what to save
        integrator  = steep  ; Algorithm (steep = steepest descent minimization)
        emtol       = 1000.0 ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
        emstep      = 0.01   ; Minimization step size
        nstenergy   = 500    ; save energies every 1.0 ps, so we can observe if we are successful
        nstxout-compressed       = 500    ; for writing coords (x) 
        nsteps      = -1     ; run as long as we need
        ; Settings that make sure we run with parameters in harmony with the selected force-field
        constraints             = none          ; no constraints
        rcoulomb                = 10            ; short-range electrostatic cutoff (in nm)
        rvdw                    = 10            ; short-range van der Waals cutoff (in nm)
        coulombtype             = Cut-Off       ; Cutoff electrostatics with large radii
        rlist                   = 20            

The main parameter here is `integrator`, which is set to `steep` or `steepest descent <https://manual.gromacs.org/current/reference-manual/algorithms/energy-minimization.html>`_ energy minimization algorithm.
Also note that we set all the cut-offs to an arbitrary large value, which is only allowed since our simulation box is large.
This is to minimize the effect of switching the interaction potential off near the cut-off radius.