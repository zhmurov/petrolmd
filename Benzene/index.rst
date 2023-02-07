Create topology and coordinates for benzene molecule
====================================================

For the molecular system that contains a mixture of small molecules it is convinient to have topology and coordinate files for each component.
This way the system coordinates can be created by a third-party software (e.g. Packmol) and the topology file will list the molecules and their quantities.
We will use an example of benzene to explain the procedure.

The coordinates
---------------

We will write basic code to syntesize the coordinates according to the atomic parameters.
Another approach would be to find the coordinates on the internet, as a stand-alone file or part of another molecular system.
The later for benzene may be, e.g. a part of phenilalanine side-chain form any protein that contains this amino-acid.
Note that the the carbon, that connects the side-chain to the backbone of the protein will have to be substituted with a hydrogen.
Also, if the system is syntesized with the correct force-field parameters, there is no need to minimize the energy (but we will do it to show the complete procedure).

Let us use CHARMM 36 force-field for the reference. If we find benzene molecule there (cgenff.rtp file), we see that all the carbons have type `CG2R61`, hydrogens are of type `HGR61`.

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

As one cas see from the comments, these atom types are found in phenilalanine and tyrosine amino-acids as well as in benzene molecule.
The equilibrium bond distances are 0.1375 nm for bond between carbons and 0.108 for carbon-hydrogen bond. The equilibrium angle is 120 degrees for both C-C-C and C-C-H angles.

So we need to creare a regular hexagonal ring with carbons in its corners and hydrogens sticking out.
Each side of the ring is a=0.1375 nm and hydrogens are farther b=0.108 nm away from the corners.

If placed in the x-y plane so that the center of hexagon is (0,0,0) and the longest diagonal along x-axis, the coordinates of the verices are (clockwise):

    .. math::

        (a,0,0), (a/2,-cos(30)*a,0), (-a/2,-cos(30)*a,0), (-a,0,0), (-a/2,cos(30)*a,0), (a/2,cos(30)*a,0)

The coordinates of the respected hydrogens are in the same layout, but (a+b) away from the center, i.e.:

    .. math::

        (a+b,0,0), ((a+b)/2,-cos(30)*(a+b),0), (-(a+b)/2,-cos(30)*(a+b),0), (-(a+b),0,0), (-(a+b)/2,cos(30)*(a+b),0), ((a+b)/2,cos(30)*(a+b),0)
