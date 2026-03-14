#!/bin/bash

# GPAW save file is too large to be included in Git repo. Therefore, to run this test, an actual GPAW calculation will be performed to obtain the calculator object.
# To install GPAW, please follow the instructions on https://gpaw.readthedocs.io/install.html

# Run GPAW to generate calculator object
mpirun python run_gpaw.py

# Run PyZFS to compute the ZFS tensor
mpirun pyzfs --wfcfmt gpaw --gpwfile nc62.gpw --ae True > zfs.out
# An equivalent command is:
# mpirun python -m pyzfs.run --wfcfmt gpaw --gpwfile nc62.gpw --ae True > zfs.out

D=`grep --color=never "D unit" zfs.xml | grep --color=never -Eoh '[+-]?[0-9]+([.][0-9]+)?'`
Dref=`grep --color=never "D unit" zfs_ref.xml | grep --color=never -Eoh '[+-]?[0-9]+([.][0-9]+)?'`

E=`grep --color=never "E unit" zfs.xml | grep --color=never -Eoh '[+-]?[0-9]+([.][0-9]+)?'`
Eref=`grep --color=never "E unit" zfs_ref.xml | grep --color=never -Eoh '[+-]?[0-9]+([.][0-9]+)?'`

echo "D = " $D
echo "Ref D = " $Dref
echo "E = " $E
echo "Ref E = " $Eref

if [ `python -c "print(int(abs($D - $Dref) < 1 and abs($E - $Eref) < 1))"` -ne 0 ]
then
    exit 0
else
    exit 1
fi

