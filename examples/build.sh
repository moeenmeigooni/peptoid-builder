#!/bin/bash

# Example build script R20 peptoid system

PYTHONDIR="/path/to/python/bin"
AMBERDIR="/path/to/amber/bin"

sequence="NMC,R20,R20,R20,R59,NH2"
replicate_ids="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
stages="6,9,12"

# 1. Prepare tleap input files for structure generation
echo "(1/10) Preparing tleap files for peptoid structure generation..."
$PYTHONDIR/python write_leap_input.py -seq $sequence -reps $replicate_ids -mode generate

# 2. Run tleap input files for structure generation
echo "(2/10) Running tleap to generate initial structure..."
IFS=',' read -ra arr <<< "$replicate_ids"
for i in "${arr[@]}"; do
    $AMBERDIR/tleap -f build_rep${i}.leap > build_rep${i}.log
done

# 3. Minimize all peptoid structures in vacuum
echo "(3/10) Minimizing structures..."
IFS=',' read -ra arr <<< "$replicate_ids"
for i in "${arr[@]}"; do
    $PYTHONDIR/python minimize.py $i >& min_rep${i}.log
done

# 4. Check for ring piercings
echo "(4/10) Checking for ring piercings..."
IFS=',' read -ra arr <<< "$replicate_ids"
for i in "${arr[@]}"; do
    $PYTHONDIR/python long_bond_checker.py -pdb peptoid_rep${i}_min.pdb -top peptoid_rep${i}.parm7 >& long_bonds.log
done

# 5. Unpierce rings, if any
echo "(5/10) Unpiercing rings, if any..."
$PYTHONDIR/python unpierce.py >& unpierce.log

# 6. Orient peptoid S-S axis along Z, determine padding distance for solvent
echo "(6/10) Orienting along S-S axis..."
vmd -dispdev text -e orient.tcl -args $replicate_ids >& orient.log

# 7. Prepare tleap input files for solvation in TIP3P
echo "(7/10) Preparing tleap files for peptoid solvation..."
$PYTHONDIR/python write_leap_input.py -seq $sequence -reps $replicate_ids -mode solvate

# 8. Run tleap to solvate all peptoid structures in TIP3P
echo "(8/10) Running tleap to solvate..."
IFS=',' read -ra arr <<< "$replicate_ids"
for i in "${arr[@]}"; do
    $AMBERDIR/tleap -f solvate_rep${i}.leap > solvate_rep${i}.log
done

# 7. Translate coords for openmm-style pbc and fix pbc dims
echo "(9/10) Fixing pbc and translating coords..."
vmd -dispdev text -e translate.tcl -args $replicate_ids >& translate.log

# 8. Re-minimize and equilibrate solvated peptoid
echo "(10/10) Equilibrating solvated peptoids..."
IFS=',' read -ra arr <<< "$replicate_ids"
for i in "${arr[@]}"; do
    IFS=',' read -ra arr2 <<< "$stages"
    for j in "${arr2[@]}"; do
        echo " >>>> Performing mineq for replicate ${i} stage ${j}"
        $PYTHONDIR/python mineq.py $i $j >& eq_rep${i}_hold${j}.log
    done
done

echo "Done"
