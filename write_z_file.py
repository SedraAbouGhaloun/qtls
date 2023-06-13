import pandas as pd
import sys


var = pd.read_csv('tmpdir/ukb22828_c'+sys.argv[1]+'_b0_v3_'+sys.argv[2]+'_'+sys.argv[3],
                  comment='#',sep='\t')

var= var.drop(['number_of_alleles','alternate_ids'], axis = 1)

var = var.rename(
        columns={'first_allele': 'allele1',
                 'alternative_alleles': 'allele2'}
    )

var = var[['rsid', 'chromosome', 'position' ,'allele1' ,'allele2']]
var.to_csv('tmpdir/ukb22828_c'+sys.argv[1]+'_b0_v3_'+sys.argv[2]+'_'+sys.argv[3]+ '.z', sep=' ', index=False)

