#!/bin/bash -l

#$ -S /bin/bash

#$ -cwd

#$ -N p3m1_msd

#$ -l h_rt=24:00:00

press=(0 50 100 200 300)

temp=(300 600 900 1200)

for i in "${press[@]}"; do

cd $i

for j in "${temp[@]}"; do

cd $j

cp ../../MSDP.py MSDP.py

python3 MSDP.py

cp MSD-overall.png ../../../../MDAnalysis/MSD/P-3m1/${i}/${j}/MSD_p3m1-overall_p${i}_t${j}.png
cp MSD-1.png ../../../../MDAnalysis/MSD/P-3m1/${i}/${j}/MSD_p3m1_1_p${i}_t${j}.png
cp MSD-2.png ../../../../MDAnalysis/MSD/P-3m1/${i}/${j}/MSD_p3m1_2_p${i}_t${j}.png
cp MSD-3.png ../../../../MDAnalysis/MSD/P-3m1/${i}/${j}/MSD_p3m1_3_p${i}_t${j}.png

cd ..

done

cd ..

done
