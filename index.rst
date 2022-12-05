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

Optionally, one can set the installation path by adding ``-DCMAKE_INSTALL_PREFIX=/your/desired/path`` to the ``cmake`` command.

3. You can optionally check the build before installing by running the tests:

    .. code-block:: shell
        
        make check -j4

4. It is convenient to set a variable for GROMACS executable. If you install into the default location:

    .. code-block:: shell

        GMX=/usr/local/gromacs/bin/gmx

Force-field and additional software:
====================================

TraPPE-UA forcefield
--------------------

Get the forcefield and make it available to GROMACS. Assuming that GROMACS in installed at ``/usr/local/gromacs``:

    .. code-block:: shell

        git clone https://github.com/wesbarnett/trappeua.git

``pdb2gmx`` utility uses enumeration for the water model selection. Hence, one need the water mode to have one of the pre-defined names, which is not the case in the ``trappeua`` implementation we just downloaded. In order to use pdb2gmx later on, we now need to rename the water model file for tip4p model and change this name in the watermodels.dat file before we copy forcefield to the GROMACS installation folder:

    .. code-block:: shell
        mv trappeua/trappeua.ff/tip4p2005.itp trappeua/trappeua.ff/tip4p.itp
        sed -i 's/tip4p2005/tip4p/' trappeua/trappeua.ff/watermodels.dat

Note, that the patched version of the forcefield is also available at:

    .. code-block:: shell
        git clone git@github.com:zhmurov/trappeua.git

Now, copy the result to GROMACS installation folder:

    .. code-block:: shell
        sudo cp -pr trappeua/trappeua.ff /usr/local/gromacs/share/gromacs/top/

CHARMM36 forcefield
-------------------

Alkanes topologies and coordinates are not standard in CHARMM forcefield. Hence, we need to use the patched version of the forcefield. The installation is similar to TraPPE-UA:

   .. code-block:: shell

      git clone git@gitlab.com:artemzhmurov/charmm36.git
      sudo cp -pr charmm36-jul2021.ff /usr/local/gromacs/share/gromacs/top/


PackMol
-------

To create the coordinates for a box of molecules, we can use Packmol software. You will need ``gfortran``, which you can install by running ``sudo apt install gfortran``. To get and install Packmol:

    .. code-block:: shell

        git clone https://github.com/m3g/packmol.git
        cd packmol
        git checkout v20.3.5
        ./configure
        make
        PACKMOL=$(pwd)/packmol



.. toctree::
   :maxdepth: 1
   :caption: Table of contents

   CreateAlkanes/index
   CreateAlkanesCHARMM/index