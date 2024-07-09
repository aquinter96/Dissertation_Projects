#!/bin/bash

module unload r
module load r/4.3.2

for datset in {1..500}
do
for sampsize in 50 60 70 80
do

sbatch -p general -N 1 --mem=32g -n 75 --cpus-per-task 1 -t 168:00:00 R CMD BATCH --vanilla --slave "--args $datset $sampsize" Model4_Simulation_Code.R

done
done