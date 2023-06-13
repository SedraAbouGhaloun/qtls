#!/bin/sh

## script to create a LD-matrix
## Sedra Abou Ghaloun 09/06/2023

## export location of files

export dir=/sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448

## get the chromosome
export chr=${1}
export lowpos=${2}
export uppos=${3}

echo "Chromosome ${chr} : Locus start ${lowpos} : Locus end ${uppos}"

## create subset BGEN file (rsids of interest)
/sc-projects/sc-proj-dh-ukb-intergenics/deps/bgenix \
-g /sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448/ukb22828_c${chr}_b0_v3.bgen \
-incl-range ${chr}:${lowpos}-${uppos} > tmpdir/ukb22828_c${chr}_b0_v3_${lowpos}_${uppos}.bgen


## create corresponding index file
/sc-projects/sc-proj-dh-ukb-intergenics/deps/bgenix \
-index -g tmpdir/filtered.${chr}.${lowpos}.${uppos}.bgen 


## run bgenix and python script to creat z-file
/sc-projects/sc-proj-dh-ukb-intergenics/deps/bgenix \
-g /sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448/ukb22828_c10_b0_v3.bgen\
-list -incl-range ${chr}:${lowpos}-${uppos} > tmpdir/ukb22828_c${chr}_b0_v3_${lowpos}_${uppos}

python write_z_file.py ${chr} ${lowpos} ${uppos}

## create master file 
python write_master_file.py ${chr} ${lowpos} ${uppos}

## run LD correlation
/sc-projects/sc-proj-dh-ukb-intergenics/deps/ldstore_v2.0_x86_64 \
--in-files tmpdir/master.z \
--write-text \
--n-threads 30 \
--read-only-bgen


## remove BGEN no longer needed
rm tmpdir/*



