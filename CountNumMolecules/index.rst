GROMACS simulations of a box of 1000 octane molecules
=====================================================

System coordinates
------------------

To create the system, we will use Packmol software.
Packmol takes coordinates of a single molecule and fills the specified volume with translated and rotated copies of the molecule.
We assume that there is a ``coord`` folder with coordinates of single molecules that can be used to create a mixture.
To create a 10nm x 10nm x 10nm box containing 1000 octane molecules with Packmol, create ``C8H18_1000_Packmol.inp`` file with the following:

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
    
        $PACKMOL < C8H18_1000_Packmol.inp

This should create a ``C8H18_1000.pdb`` file. Feel free to load it into VMD or other visualization software to have a look.

System topology
---------------

Now we have coordinates for the system, we need to create topology file.
Since we already created the topology files for all the molecules in the mixture (octane in this case), the system topology is quite simple.
We assume that these files are located in ``toppar`` folder in your working directory.
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
        $GMX solvate -cp C8H18_1000_box.gro -o C8H18_1000_solv.gro -p C8H18_1000.top

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

You will still need to run ``gmx editconf`` to specify the box size:

    .. code-block:: shell
    
        $GMX editconf -f C8H18_1000_solv.pdb -o C8H18_1000_solv.gro -box 10 10 10

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

Note that the name of water molecule in the topology is ``SOL``, which is how it is called in the forcefield by default.
The drawback of this approach is that user has to pre-compute the number of water molecules instead of relying on GROMACS to add water up to desired density.
The first advantage is that we eliminate one step in the procedure of system generation.
But more importantly, we can slightly modify the Packmol input file and create molecular system where the two components are separated.
This is done by providing different compartments to Packmol:

    .. code-block:: text

        tolerance 2.0
        filetype pdb
        output C8H18_1000_solv.pdb

        structure coord/C8H18.pdb
        number 1000 
        inside box 0. 0. 0. 50. 100. 100. 
        end structure

        structure coord/tip3p.pdb
        number 5000 
        inside box 50. 0. 0. 100. 100. 100. 
        end structure

This way, the simulation box will be split in half with the left side filled with octane and right side filled with water.
Note that Packmol provides an extensive number of options to specify the geometry of the system.
See `Packmol users guide <https://m3g.github.io/packmol/userguide.shtml>`_ for details.

Preparing the system
--------------------

Energy minimization
^^^^^^^^^^^^^^^^^^^

    The configuration file for energy minimization follows.

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
        constraints             = h-bonds   ; bonds involving H are constrained
        rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
        rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
        vdw-modifier            = Force-switch ;  specific CHARMM
        rvdw_switch             = 1.0       ;
        DispCorr                = no        ; account for cut-off vdW scheme -
        ;in case of CHARMM DispCorr = EnerPres only for monolayers
        coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
        fourierspacing          = 0.15     ; grid spacing for FFT

This file is very similar to the one we used for vacuum minimization with only non-bonded parameters different.
The reason is that the box is now not infinitely large and cut-offs should be adjusted to fit this size.
The parameters used here are recommended parameters for CHARMM force-fields, i.e. the force-field was parametrize with these non-bonded parameters in mind.

    .. code-block:: shell

        $GMX grompp -f em.mdp -c C8H18_1000_solv.gro -p C8H18_1000.top -o em.tpr
        $GMX mdrun -deffnm em

Equilibration
^^^^^^^^^^^^^

The equilibration is usually done in two steps.
First run of equilibration is with constant volume (NVT ensemble).
This is so that the barostat will not pick up large deviation in the pressure that may occur because of bad placement of the molecules.
The configuration file for the NVT equilibration in CHARMM force-field is:

    .. code-block:: text
        
        title                   = NVT equilibration 

        ; Parameters describing what to do, when to stop and what to save
        integrator              = md        ; leap-frog integrator
        dt                      = 0.002     ; 2 fs
        nsteps                  = 50000     ; 2 * 50000 = 100 ps
        nstenergy               = 500       ; save energy and temperature every 1.0 ps
        nstxout-compressed      = 5000    ; for writing coords (x) 

        ; periodic boundary condition
        pbc                     = xyz       ;

        ; Keep system temperature fluctuating physically correct
        tcoupl                  = V-rescale           ; modified Berendsen thermostat
        tc-grps                 = system   ; coupling groups 
        tau_t                   = 0.1      ; time constant, in ps
        ref_t                   = 300      ; reference temperature, one for each group, in K

        ; Pressure coupling is off
        pcoupl                  = no

        ; Velocity generation
        gen_vel                 = yes                 ; assign velocities from Maxwell distribution
        gen_temp                = 300                 ; temperature for Maxwell distribution

        ; Settings that make sure we run with parameters in harmony with the selected force-field
        constraints             = h-bonds   ; bonds involving H are constrained
        rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
        rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
        vdw-modifier            = Force-switch ;  specific CHARMM
        rvdw_switch             = 1.0       ;
        DispCorr                = no        ; account for cut-off vdW scheme -
        ;in case of CHARMM DispCorr = EnerPres only for monolayers
        coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
        fourierspacing          = 0.15     ; grid spacing for FFT

Here, we use ``md`` integrator, which is a leap-frog scheme.
Note that now we save the coordinates every 5000 steps in compressed format, so that we can monitor the progress (``nstxout-compressed`` parameter).
We employ the temperature control with `velocity rescaling algorithm <https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html#mdp-tcoupl>`_.
The pressure coupling is off on this stage.
The initial velocities are generated at this stage, based on the temperature of 300K.
Non-bonded parameters remain the same as they should for the rest of the simulations.
To start the NVT equilibration stage, create portable simulation file (``.tpr``) with ``gmx grompp`` and start the run with ``gmx mdrun``.


    .. code-block:: shell

        $GMX grompp -f nvt.mdp -c em.gro -p C8H18_1000.top -o nvt.tpr
        $GMX mdrun -deffnm nvt

The second equilibration stage is done under constant pressure conditions (NPT ensemble).
GROMACS will be adjusting the size of the simulation box to reach the target value of pressure.
The volume of the box can change quite drastically at this stage if the initial box is overfilled or underfilled.
We aim to reach more or less conserved volume at the end of this run as an indicator that the system is well equilibrated and ready for production run.

    .. code-block:: text

        title                   = NPT equilibration 

        ; Parameters describing what to do, when to stop and what to save
        integrator              = md        ; leap-frog integrator
        dt                      = 0.002     ; 2 fs
        nsteps                  = 50000     ; 2 * 50000 = 100 ps
        nstenergy               = 500       ; save energy and temperature every 1.0 ps
        nstxout-compressed      = 5000    ; for writing coords (x) 

        ; periodic boundary condition
        pbc                     = xyz       ;

        continuation            = yes      

        ; Pressure coupling is on
        pcoupl                  = C-rescale             ; Pressure coupling on in NPT
        pcoupltype              = isotropic             ; uniform scaling of box vectors
        tau_p                   = 1.0                   ; time constant, in ps
        ref_p                   = 1.0                   ; reference pressure, in bar
        compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
        refcoord_scaling        = com

        ; Keep system temperature fluctuating physically correct
        tcoupl                  = V-rescale           ; modified Berendsen thermostat
        tc-grps                 = system   ; coupling groups 
        tau_t                   = 0.1      ; time constant, in ps
        ref_t                   = 300      ; reference temperature, one for each group, in K

        ; Settings that make sure we run with parameters in harmony with the selected force-field
        constraints             = h-bonds   ; bonds involving H are constrained
        rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
        rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
        vdw-modifier            = Force-switch ;  specific CHARMM
        rvdw_switch             = 1.0       ;
        DispCorr                = no        ; account for cut-off vdW scheme -
        ;in case of CHARMM DispCorr = EnerPres only for monolayers
        coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
        fourierspacing          = 0.15     ; grid spacing for FFT

Here we have isotropic pressure coupling enabled with exponential relaxation pressure coupling scheme.
More on pressure coupling algorithms that are supported by GROMACS `here <https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html#mdp-pcoupl>`_.
Now, configure and run the simulations.

    .. code-block:: shell

        $GMX grompp -f npt.mdp -c nvt.gro -p C8H18_1000.top -o npt.tpr
        $GMX mdrun -deffnm npt


Production simulations
----------------------

Production simulations can be berformed with isotropic or anisotropic pressure coupling scheme.
Configuration file for isotropic pressure coupling is very similar to the one for the NPT equilibration.
All we need to do is to adjust number of steps and how often we want the output to be saved.

    .. code-block:: text

        title                   = Equilibrium simulations

        ; Parameters describing what to do, when to stop and what to save
        integrator              = md        ; leap-frog integrator
        dt                      = 0.002     ; 2 fs
        nsteps                  = 5000000     ; 

        ; periodic boundary condition
        pbc                     = xyz       ;

        continuation            = yes      

        ; Output control - output frequency in steps
        ; Output frequency for  output trajctory file ,trr
        nstxout                  = 0       ; for writing coords (x) 
        nstvout                  = 0       ; for writing velocities (v) 
        nstfout                  = 0       ; for writing forces (f) 
        ; Output frequency for energies to log file and energy file
        nstlog                   = 1000    ; for writing energies to log file
        nstenergy                = 500     ; for writing energies to edr file 
        ; Output frequency and precision for .xtc file
        nstxout-compressed       = 5000    ; for writing coords (x) 

        ; Pressure coupling is on
        pcoupl                  = C-rescale             ; Pressure coupling on in NPT
        pcoupltype              = isotropic             ; uniform scaling of box vectors
        tau_p                   = 5.0                   ; time constant, in ps
        ref_p                   = 1.0                   ; reference pressure, in bar
        compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
        refcoord_scaling        = com

        ; Keep system temperature fluctuating physically correct
        tcoupl                  = V-rescale           ; modified Berendsen thermostat
        tc-grps                 = system ; 
        tau_t                   = 0.1    ; time constant, in ps
        ref_t                   = 300    ; reference temperature, one for each group, in K

        ; Settings that make sure we run with parameters in harmony with the selected force-field
        constraints             = h-bonds   ; bonds involving H are constrained
        rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
        rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
        vdw-modifier            = Force-switch ;  specific CHARMM
        rvdw_switch             = 1.0       ;
        DispCorr                = no        ; account for cut-off vdW scheme -
        ;in case of CHARMM DispCorr = EnerPres only for monolayers
        coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
        fourierspacing          = 0.15     ; grid spacing for FFT

For anisotropic pressure coupling, we will have to use Parrinello-Rahman barostat and set the target pressure for each component individually.
The box will shrink in different dimentions differently, depending on the pressure along the respective component.
This type of pressure coupling is usefull, if the system itself is anisotropic (e.g. when there is a flat layer of compound or clear separation between mixed component along one of the axis).
To set up anisotropic pressure coupling, use the following block in ``.mdp`` file:

    .. code-block:: text

        ; Pressure coupling is on
        pcoupl                  = Parrinello-Rahman                         ; Pressure coupling on in NPT
        pcoupltype              = anisotropic                               ; non-uniform scaling of box vectors
        ref_p                   = 1.0    1.0    1.0    0.0    0.0    0.0    ; reference pressure, in bar. No shear, off-diagonal elements are zero
        tau_p                   = 50.0                                       ; time constant, in ps
        compressibility         = 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 ; isothermal compressibility of water, bar^-1
        refcoord_scaling        = com

Save the ``.mdp`` file, create ``.tpr`` with ``gmx grompp`` and run the simulations with ``gmx mdrun``:

    .. code-block:: shell

        $GMX grompp -f md.mdp -c npt.gro -p C8H18_1000.top -o md.tpr
        $GMX mdrun -deffnm md

