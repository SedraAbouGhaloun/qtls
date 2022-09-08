import numpy as np 
import os
import pandas as pd
from scipy.sparse import coo_matrix
import torch     
import multiprocessing as mp  

def preprocess_subset(subset, prev_rs_id_list, prev_variant_list, prev_mtoid_list, prev_mtid_list, prev_study_list):
    
    """ maps every new molecular trait object id, variant,rs_id or study to a new index when it does not exist in 
     prev_mtid_list, prev_variant_list, prev_study_list, and ignores already exixting ones. returns np array of 
     the modified subset and np array of the p values """
    
    # update lists with new values and creat dictionary 
    
    prev_mtoid_list += list(set(subset['molecular_trait_object_id']) - set(prev_mtoid_list))
    mtoid_dict =  {num:i for i,num in enumerate(prev_mtoid_list)}
    
    prev_mtid_list += list(set(subset['molecular_trait_id']) - set(prev_mtid_list))
    mtid_dict =  {num:i for i,num in enumerate(prev_mtid_list)}
    
    prev_study_list += list(set(subset['study']) - set(prev_study_list))
    study_dict =  {num:i for i,num in enumerate(prev_study_list)}
    
    prev_variant_list += list(set(subset['variant']) - set(prev_variant_list))
    variant_dict =  {num:i for i,num in enumerate(prev_variant_list)}
    
    prev_rs_id_list += list(set(subset['rsid']) - set(prev_rs_id_list))
    rs_dict =  {num:i for i,num in enumerate(prev_rs_id_list)}

    # select only relevant columns
    selected_data = subset[['molecular_trait_object_id', 'molecular_trait_id', 'variant', 'pvalue','rsid','study']]

    # do the mapping
    molecular_trait_id_list = (selected_data["molecular_trait_id"].apply(lambda x: mtid_dict[x])).tolist()
    molecular_trait_object_id_list  = (selected_data["molecular_trait_object_id"].apply(lambda x: mtoid_dict[x])).tolist()
    variant_list = (selected_data["variant"].apply(lambda x: variant_dict[x])).tolist()
    study_list = (selected_data["study"].apply(lambda x: study_dict[x])).tolist()
    rsid_list = (selected_data["rsid"].apply(lambda x: rs_dict[x])).tolist()


    return rsid_list, variant_list, molecular_trait_object_id_list, molecular_trait_id_list, study_list, subset.pvalue.tolist()



def filter_subset(subset, gwas_data, p_threshold):
    
    """filters a subset based on their presence on gwas_data and their p value, returns filtered dataframe"""
    
    merged = pd.merge(subset, gwas_data, left_on='rsid', right_on='ID',indicator = True)
    return(merged[merged['pvalue'] < p_threshold])



# Read the gwas data 

gwas_path = '/dh-projects/uk_bb_intergenics/analysis/development/upmeiejv/results/gwas-general/pvar_externded_ids.tsv'
gwas_data = pd.read_csv(gwas_path,  sep= '\t',usecols= ['ID'])

# hier we save the unique molecular_trait_object_ids, molecular_trait_ids, variants and studies. 
prev_rs_id_list, prev_mtoid_list, prev_mtid_list, prev_variant_list, prev_study_list = [],[],[],[],[]

# hier we save the coordinates of the p values 
rs_id_list, molecular_trait_object_id_list, molecular_trait_id_list, variant_list, study_list, pvalue_list = [],[],[],[],[],[]

# global paramters

p_threshold = 0.005

# list all studies 

dir_path  = '.'
studies = []
for file in os.listdir(dir_path):
    if file.endswith('.tsv'):
        studies.append(file)
   
        
for study in studies:
    df =  pd.read_csv(study, sep='\t',usecols =['molecular_trait_object_id', 'molecular_trait_id', 'variant', 'pvalue','rsid'])
    print('reading done')
    filtered = filter_subset(df, gwas_data, p_threshold) 
    print('filter done')
    filtered['study'] = study.split('.')[:-2][0]
    new_rsid_list, new_variant_list, new_molecular_trait_object_id_list, new_molecular_trait_id_list,new_study_list, new_pvalue = preprocess_subset(filtered, prev_rs_id_list, prev_variant_list, prev_mtoid_list, prev_mtid_list, prev_study_list)   
    print('preprocess done')
    molecular_trait_object_id_list += new_molecular_trait_object_id_list
    molecular_trait_id_list += new_molecular_trait_id_list
    variant_list += new_variant_list
    study_list += new_study_list
    pvalue_list += new_pvalue
    
    

with open("molecular_trait_object_id.csv", "w") as f:
    f.write('\t'.join(map(str, molecular_trait_object_id_list)))

with open("molecular_trait_id.csv", "w") as f:
    f.write('\t'.join(map(str, molecular_trait_id_list)))
    
with open("variant.csv", "w") as f:
    f.write('\t'.join(map(str, variant_list)))
         
with open("rsid.csv", "w") as f:
    f.write('\t'.join(map(str, rs_id_list)))


with open("study.csv", "w") as f:
    f.write('\t'.join(map(str, study_list)))    

with open("pvalue.csv", "w") as f:
    f.write('\t'.join(map(str, pvalue_list)))
    
with open("molecular_trait_object_id_idx.csv", "w") as f:
    f.write('\t'.join(map(str,prev_mtoid_list)))
    
with open("molecular_trait_id_idx.csv", "w") as f:
    f.write('\t'.join(map(str,prev_mtid_list)))    
    
with open("variant_idx.csv", "w") as f:
    f.write('\t'.join(map(str,prev_variant_list)))

    
with open("rsid_idx.csv", "w") as f:
    f.write('\t'.join(map(str,prev_rs_id_list )))
    
with open("study_idx.csv", "w") as f:
    f.write('\t'.join(map(str,prev_study_list)))
    


