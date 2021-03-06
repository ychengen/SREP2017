#!/bin/bash

module load nco

## declare an array storing ESM names
declare -a arr=("BNU-ESM" "HadGEM2-CC" "HadGEM2-ES" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM" "MIROC-ESM-CHEM")

for i in "${arr[@]}"
   do
      if [ "$i" = "HadGEM2-CC" -o "$i" = "HadGEM2-ES" ]; then
        sed -s "s/modelnm/${i}/g"  template_PI2.sh         > run_avg_${i}_PI.sh
      else
        sed -s "s/modelnm/${i}/g"  template_PI1.sh         > run_avg_${i}_PI.sh
      fi
      chmod 755 ./run_avg_${i}_PI.sh
      ./run_avg_${i}_PI.sh
      rm -f ./run_avg_${i}_PI.sh
done

wait

