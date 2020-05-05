#!/bin/bash -l

#$ -S /bin/bash$!/bin.bash

temp=(1500 1800)
press=(0 50 100 200 300)

for p in "${press[@]}"; do

cd $p

for t in "${temp[@]}"; do

mkdir $t

cd $t

cp ../../POTCAR POTCAR
cp ../../KPOINTS KPOINTS
cp ../init_${p}.POSCAR.vasp POSCAR
cp ../../INCAR INCAR
cp ../../run.bash run.bash

sed "/TEBEG/s/.*/TEBEG = ${t}/g" INCAR > INCtemp
rm INCAR
sed "/TEEND/s/.*/TEEND = ${t}/g" INCtemp > INCAR
rm INCtemp
sed "/#$ -N p3m1/s/.*/#$ -N p3m1p${p}t${t}/g" run.bash > runtemp.bash
mv runtemp.bash run.bash

cd ..

done

cd ..

done
