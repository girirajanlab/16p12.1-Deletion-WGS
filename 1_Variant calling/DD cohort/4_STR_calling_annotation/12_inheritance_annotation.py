import pandas as pd

# Add inheritance annotations for STRs

# Input and output files
STRS="/path/to/strs.csv" # Use the output of script 11_annotate_loeuf.py
FAM="/path/to/sample/FAM/file.fam"

OUTPUT_FILE="/path/to/output.csv"

# Load files
strs = pd.read_csv(STRS, index_col = 0)
fam = pd.read_csv(FAM, sep = '\t')

def add_inheritance(row, mother = False, father = False):
    # Get parent information
    sample = row['Sample']
    
    mother_code = fam[fam.IID==sample]['Mother'].to_string(index = False, header = False).strip()
    father_code = fam[fam.IID==sample]['Father'].to_string(index = False, header = False).strip()
    
    # Get carrier status of parents
    if mother_code == '0':
        mother_return = '0'
    elif fam[fam.IID==mother_code]['Carrier'].to_string(index = False, header = False).strip()=='unknown':
        mother_return = 'M'
    elif int(fam[fam.IID==mother_code]['Carrier'])==1:
        mother_return = 'MC'
    elif int(fam[fam.IID==mother_code]['Carrier'])==0:
        mother_return = 'MNC'
    else:
        mother_return = 'M'

    if father_code == '0':
        father_return = '0'
    elif fam[fam.IID==father_code]['Carrier'].to_string(index = False, header = False).strip()=='unknown':
        father_return = 'F'
    elif int(fam[fam.IID==father_code]['Carrier'])==1:
        father_return = 'FC'
    elif int(fam[fam.IID==father_code]['Carrier'])==0:
        father_return = 'FNC'
    else:
        father_return = 'F'
        
    # Check if parent has call
    vid = row['variant_id']
    vid_samps = strs[strs.variant_id==vid]['Sample'].to_list()
    
    if mother_code in vid_samps:
        mother = True
    if father_code in vid_samps:
        father = True
    
    if mother and father:
        return '.'
    if mother:
        return mother_return
    if father:
        return father_return    
    return '.'
    
strs['inheritance'] = strs.apply(add_inheritance, axis = 1)

strs.to_csv(OUTPUT_FILE, index = False)