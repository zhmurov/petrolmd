Computing surface tension for hydrocarbons
==========================================

Calculating pressure from MD
----------------------------

The surface tension can be computed using the following relation [AlejandreJCP95]_:

    .. math::

        \gamma=\frac{L_{X}}{2}\langle P_{XX}-\frac{1}{2}\left(P_{YY}+P_{ZZ}\right)\rangle




.. [AlejandreJCP95] Alejandre, J., Tildesley, D. J., & Chapela, G. A. (1995). `Molecular dynamics simulation of the orthobaric densities and surface tension of water. <https://aip.scitation.org/doi/pdf/10.1063/1.469505>`_ J. Chem. Phys., 102(11), 4574--4583.