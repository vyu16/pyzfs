cd o2_qbox_xml
pyzfs --wfcfmt qbox --filename o2.xml > zfs.out
cd ..

cd o2_qe_hdf5
pyzfs --wfcfmt qeh5 --memory high > zfs.out
cd ..

cd c3h6_gpaw
python3 run_gpaw.py
pyzfs --wfcfmt gpaw --gpwfile c3h6.gpw > zfs.out
cd ..

cd c3h6_gpaw_ae
python3 run_gpaw.py
pyzfs --wfcfmt gpaw --gpwfile c3h6.gpw --ae True > zfs.out
cd ..
