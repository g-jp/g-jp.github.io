
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter

# This collection was created for data analysis regarding mass spectrometry and the functions described are to be used along with 
# with the MassTRIX data sheet not comparing isotopes. 

# All functions' syntaxes follow the same principle: function(KEGG_cids_list, File_1, File_2), in which KEGG_cids_list is a list with
# all the KEGG_cids from intended compound, perhaps obtained from KEGG_mapper tool for metabolism comparison, File_1 is 
# the name of the file refering to the first organism to be compared (usually the reference one) and File_2 is the name of 
# the second file refering to the second organism to be compared

# For specific signal attribution, make sure the data sheet obtained from MassTRIX has most 1 ppm or less since a sign for the
# ionization pattern for a molecule may be attributed to two peaks incorrectly.

# ----------------------------------------------------------------------------------------------------

def plot_cids(x, y, z): #Creates a plot for a list of KEGG_cids showing comparative values of normalized peak height for metabolites in both BY and KO types
# ---------------------------------------------
    Type1 = pd.read_csv(str(y), sep = '\t') 
    Type2 = pd.read_csv(str(z), sep = '\t')
    T1allmet_one = Type1.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid')
    T1allmet = Type1.set_index('KEGG_cid')
    T2allmet_one  = Type2.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid')
    T2allmet = Type2.set_index('KEGG_cid')
    cid_listBY = []
    for n in Type1.KEGG_cid:
        cidBY = n
        cid_listBY.append(cidBY)
    cid_listKO = []
    for n in Type2.KEGG_cid:
        cidKO = n
        cid_listKO.append(cidKO)
# ---------------------------------------------

    ord = []
    for n in x:
        if n in cid_listBY and n in cid_listKO:
            val = []
            heighta = T2allmet.peak_height[str(n)].sum()
            heightb = T1allmet.peak_height[str(n)].sum()
            enkBY = T1allmet.peak_height.HMDB01045.sum()
            enkKO = T2allmet.peak_height.HMDB01045.sum()
            qnt_a = heighta/enkKO
            qnt_b = heightb/enkBY
            val.append(qnt_a)
            val.append(qnt_b)
            for c in val:
                ord.append(c)
        else:
            if n in cid_listBY:
                print(f'{n} appears only in {y} type')  
            elif n in cid_listKO:
                print(f'{n} appears only in {z} type')
            else:
                print(f'{n} does not appear as a signal')
    abc = []
    for n in x:
        if n in cid_listBY and n in cid_listKO:
            name = []
            a = T2allmet_one.KEGG_name[str(n)].split(';')[0] 
            b = T1allmet_one.KEGG_name[str(n)].split(';')[0] 
            if a[-1] == ')':
                a = a[:-9] + '\nT2'
                b = b[:-9] + '\nT1'
            else:
                a = a + '\nT2'
                b = b + '\nT1'    
            name.append(a)
            name.append(b)
            for c in name:
                abc.append(c)
        else:
            continue
    with sns.axes_style('whitegrid'):
        plt.subplots(figsize=(6,6))
        c = sns.barplot(ord,abc)
    return c

def df_cids(x, y, z): #Creates a DataFrame for metabolites in both BY and KO types, comparing the normalized peak values for each metabolite in a list
# ---------------------------------------------
    Type1 = pd.read_csv(str(y), sep = '\t') 
    Type2 = pd.read_csv(str(z), sep = '\t')
    T1allmet_one = Type1.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid')
    T1allmet = Type1.set_index('KEGG_cid')
    T2allmet_one  = Type2.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid')
    T2allmet = Type2.set_index('KEGG_cid')
    cid_listBY = []
    for n in Type1.KEGG_cid:
        cidBY = n
        cid_listBY.append(cidBY)
    cid_listKO = []
    for n in Type2.KEGG_cid:
        cidKO = n
        cid_listKO.append(cidKO)
# ---------------------------------------------

    ordKO = []
    ordBY = []
    for n in x:
        if n in cid_listBY and n in cid_listKO:
            val = []
            heighta = T2allmet.peak_height[str(n)].sum()
            heightb = T1allmet.peak_height[str(n)].sum()
            enkBY = T1allmet.peak_height.HMDB01045.sum()
            enkKO = T2allmet.peak_height.HMDB01045.sum()
            qnt_a = heighta/enkKO
            qnt_b = heightb/enkBY
            val.append(qnt_a)
            val.append(qnt_b)
            ordKO.append(val[0])
            ordBY.append(val[1])
        else:
            continue
    abcKO = []
    abcBY = []
    for n in x:
        if n in cid_listBY and n in cid_listKO:
            name = []
            a = T2allmet_one.KEGG_name[str(n)].split(';')[0] 
            b = T1allmet_one.KEGG_name[str(n)].split(';')[0] 
            if a[-1] == ')':
                a = a[:-9] + '\nKO'
                b = b[:-9] + '\nBY'
            else:
                a = a + '\nKO'
                b = b + '\nBY'    
            name.append(a)
            name.append(b)
            abcKO.append(name[0][:-3])
            abcBY.append(name[1][:-3])                
        else:
            continue  

    dictKO = dict(zip(abcKO, ordKO))
    dictBY = dict(zip(abcBY, ordBY))
    dict_list = [dictKO,dictBY]
    
    c = pd.DataFrame(dict_list, index=['T2','T1']).transpose()
    return c

def quantify(x, y, z): #Calculate the normalized value for an specific metabolite based on its KEGG_cid

# ---------------------------------------------
    Type1 = pd.read_csv(str(y), sep = '\t') 
    Type2 = pd.read_csv(str(z), sep = '\t')
    T1allmet_one = Type1.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid')
    T1allmet = Type1.set_index('KEGG_cid')
    T2allmet_one  = Type2.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid')
    T2allmet = Type2.set_index('KEGG_cid')
    cid_listBY = []
    for n in Type1.KEGG_cid:
        cidBY = n
        cid_listBY.append(cidBY)
    cid_listKO = []
    for n in Type2.KEGG_cid:
        cidKO = n
        cid_listKO.append(cidKO)
# ---------------------------------------------

    enkBY = T1allmet.peak_height.HMDB01045.sum()
    enkKO = T2allmet.peak_height.HMDB01045.sum()
    if x in cid_listBY and x in cid_listKO:
        heighta = T2allmet.peak_height[str(x)].sum()
        heightb = T1allmet.peak_height[str(x)].sum()
        qnt_a = heighta/enkKO
        qnt_b = heightb/enkBY
        name = T2allmet_one.KEGG_name[str(x)].split(';')[0]
        if name[-1] == ')':
            name = name[:-9] 
        else:
            name = name
        return print(f'The normalization for {name} is {qnt_b} for T1 and {qnt_a} for T2')
    else:
        if x in cid_listBY:
            b = T1allmet.peak_height[str(x)].sum()
            qnt_b = b/enkBY
            name = T1allmet_one.KEGG_name[str(x)].split(';')[0]
            return print(f'The normalization for {name} is {qnt_b} in T1 only')
        elif x in cid_listKO :
            a = T2allmet.peak_height[str(x)].sum()
            qnt_a = a/enkKO
            name = T2allmet_one.KEGG_name[str(x)].split(';')[0]
            return print(f'The normalization for {name} is {qnt_a} in T2 only') 
            