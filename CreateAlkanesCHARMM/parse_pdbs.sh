#!/bin/bash

for filename in *.pdb; do
    # Replace the atom types with CHARMM naming
    sed -i 's/ CH2  /CC32A/' ${filename}
    sed -i 's/ CH3  /CC33A/' ${filename}
done