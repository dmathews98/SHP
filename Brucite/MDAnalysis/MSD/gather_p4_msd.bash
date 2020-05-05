#!/bin/bash -l

#$ -S /bin/bash

#$ -cwd

#$ -N p4_msd

#$ -l h_rt=24:00:00

press=(200 300 400)

temp=(300 600 900 1200 1500 1800)

cd P4_32_12

for i in "${press[@]}"; do

cd $i

for j in "${temp[@]}"; do

cd $j

cp p4_msd_H_p${i}_t${j}.txt ../Hgraphing/p4_msd_H_p${i}_t${j}.txt

cd ..

done

cd Hgraphing

python3 graph_combined_H.py

cd ..

cd ..

done
