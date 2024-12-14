import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as ml
import pandas as pd
import h5py
import illustris_python as il
import os
import sys
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

#Gas cut parameters (see line 55)
n1 = [] 
n3 = [] 
n3s = [] #file name string for n3

#Halo Cut choice parameters
n2 = []
n4 = []
n4s= [] #file naming string for n4 

select_files = [0]
#TNG Parameter fields 
fSH = ['SubhaloMass','SubhaloMassType', 'SubhaloGrNr']
fH = ['GroupMassType','GroupMass', 'GroupLenType','GroupFirstSub']
run = '2160'
basePath = 'TNGruns_directory/Name_of_folder_with_run'+run+'/output_folder'
res = 'HR' #resolution, low, medium, or high (based on no. of particles in the TNG run)
met_str = [] #metallicity cut file naming string
#Pick seeding percentage and string (i.e, 0.01 and '01')
seeding_percentage = 0.01 #fraction of seeded black holes at each snapshot
seeding_str = '01' #fraction of seeded black holes in the form of a file naming string 
BHmass= 10**4 #set a BH mass cut 
table_range = np.arange(0,99,1).tolist() #start at 0 for MR/HR, start at 1 for LR #this is the first index of the table to read for a subhalo/halo ID table (rows) assorted by snapshot (columns) for merger trees constrained by the gas and halo cuts #the table is transposed in line 75
V = (35/0.6774)**3 #TNG Volume

#TNG Scaling relation with scatter: function & parameters from Evans et al. 2024 
alpha = 8.46
beta = 1.05
epsilon0 = 0.34
Nparams = 3
Ngal = 35
alphaTNG = 6.68
betaTNG = 0.79
eps0_TNG = 0.1354795985658518
def ms_mBH(Mbulge2):
    mBH2 = 10**(alphaTNG + betaTNG*(np.log10(Mbulge2/10**9)+np.random.normal(0.0, eps0_TNG,len(Mbulge2))))
    return mBH2

q = int(sys.argv[1]) 
met_str_ind = met_str[q]

print(q)
for k in []: #Pick parameter cuts
    files = ['mg'+str(n1[k])+str(n3s[k])+'_mH'+str(n2[k])+str(n4s[k])+'_zpt'+met_str_ind,
           'mg'+str(n1[k])+str(n3s[k])+'_mH'+str(n2[k])+str(n4s[k])+'_zpt'+met_str_ind+'sfr',
           'mg'+str(n1[k])+str(n3s[k])+'_mSH'+str(n2[k])+str(n4s[k])+'_zpt'+met_str_ind,
           'mg'+str(n1[k])+str(n3s[k])+'_mSH'+str(n2[k])+str(n4s[k])+'_zpt'+met_str_ind+'sfr',
           'mg'+str(n1[k])+str(n3s[k])+'_mcSH'+str(n2[k])+str(n4s[k])+'_zpt'+met_str_ind,
           'mg'+str(n1[k])+str(n3s[k])+'_mcSH'+str(n2[k])+str(n4s[k])+'_zpt'+met_str_ind+'sfr',
           'TNGcut_mH5e10msolperh']
    for m in select_files:
        print(files[m])
        #No. of total BHs in Model
        fm = open('number_densities_output_folder/nBHs_'+files[m]+'noRNG_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.txt','a+')
        #No. of BHs in stochastic seeding model using RNG
        fmr = open('number_densities_output_folder/nBHs_'+files[m]+'_fpt'+seeding_str+'_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.txt','a+')
        if os.path.getsize('/n/holystore01/LABS/hernquist_lab/Users/alawrence/analislawrence/python_scripts/05_27_21/'+res+'_broad_cut_run_lenientcuts/new/desc_' + files[m] + res+'.csv')>0:
            org_data = np.array(np.loadtxt('/n/holystore01/LABS/hernquist_lab/Users/alawrence/analislawrence/python_scripts/05_27_21/'+res+'_broad_cut_run_lenientcuts/new/desc_' + files[m] + res+'.csv',delimiter=' ', usecols = None))
            
        else:
            continue
        
        org_data_T = np.array(org_data).transpose()
        

        nm = ([])
        nm_rand = ([])
        BHs_already_seeded = ([])
        BHs_already_seeded_rand = ([])
        get_desc0 = np.array(org_data_T[0])
        find_new_seeds_arr = [[],[]]
        np.random.seed()
        random_SHs0 = np.random.randint(0, len(get_desc0[~np.isnan(get_desc0)]),np.round(int(seeding_percentage*len(get_desc0[~np.isnan(get_desc0)])),0))
        random_SHs_list = random_SHs0.tolist()
    
        for j in table_range:        
            
            if m in []: #indices for halo descendant files only
                subhalos = il.groupcat.loadSubhalos(basePath, j+1, fields = fSH)
                SH_Gr_Nr = np.array(subhalos['SubhaloGrNr']).astype('int32')
                halos = il.groupcat.loadHalos(basePath, j+1, fields = fH)
                head = il.groupcat.loadHeader(basePath,j+1)
                h = head['HubbleParam']
                zs = head['Redshift']
                halo_mass = np.array(halos['GroupMass'])*1e10/h
                
                if (j > table_range[0]):
                    psubhalos = il.groupcat.loadSubhalos(basePath, j, fields = fSH)
                    pSH_Gr_Nr = np.array(psubhalos['SubhaloGrNr']).astype('int32')
                    
                get_desc = np.array(org_data_T[j]).astype('int') 
                get_desc_prev = np.array(org_data_T[j-1]).astype('int') 
                masked_desc = SH_Gr_Nr[(get_desc[(get_desc>=0)&(~np.isnan(get_desc))])]
                
                if j==table_range[0]:
                    masked_desc_prev = []
                else:                     
                    masked_desc_prev0= get_desc_prev[(get_desc>=0)&(~np.isnan(get_desc))]
                    masked_desc_prev = (pSH_Gr_Nr[(masked_desc_prev0[(masked_desc_prev0>=0)&(~np.isnan(masked_desc_prev0))])]).tolist()
                    x_ind1 = np.where((masked_desc_prev0<0)&(~np.isnan(masked_desc_prev0)))[0]
                    if len(x_ind1) > 0:
                        for q in x_ind1.tolist():
                            masked_desc_prev.insert(q,-1000)
                
                get_born_desc0 = get_desc[~np.isnan(get_desc)]
                get_born_desc = (SH_Gr_Nr[(get_born_desc0[(get_born_desc0>=0)])]).tolist()
                x_ind2 = np.where(get_born_desc0<0)[0]
                if len(x_ind2) > 0:
                    for w in x_ind2.tolist():
                        get_born_desc.insert(w,-1000)
                
                if j==table_range[0]:
                    get_born_desc_prev = []
                else:
                    get_born_desc_prev0= get_desc_prev[~np.isnan(get_desc)]
                    get_born_desc_prev = (pSH_Gr_Nr[(get_born_desc_prev0[(get_born_desc_prev0>=0)])]).tolist()
                    x_ind3 = np.where(get_born_desc_prev0<0)[0]
                    if len(x_ind3) > 0:
                        for z in x_ind3.tolist():
                            get_born_desc_prev.insert(z,-1000)
                
                np.random.seed()
                random_SHs = np.random.randint(0, len(get_born_desc),np.round(int(seeding_percentage*len(get_born_desc)),0))
                random_SHs_list.extend(random_SHs)
                rs = np.array(np.unique(random_SHs_list), dtype = int)
                get_desc_rand = np.array(get_born_desc)[(rs)]
                if j == table_range[0]:
                    get_desc_prev_rand = []
                    gdr = np.array([])
                else: 
                    get_desc_prev_rand = np.array(get_born_desc_prev)[(rs)]
                    gdr = get_desc_rand[(get_desc_rand>=0)]
                if j ==table_range[0]:
                    pdr = []
                else:
                    pdr=np.array(get_desc_prev_rand[np.where(get_desc_prev_rand>=0)])
                
                if len(masked_desc) == 0:
                    print("none")
                    nm.append(0)
                    BHs_already_seeded.append(0)
                    nm_rand.append(0)
                    BHs_already_seeded_rand.append(0)
                    #save zero mass densities
                    f_rhom = open('mass_density_directory/rhom_scalingrelation_mBH_ms_'+files[m]+'_fpt'+seeding_str+'_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.csv','a+')
                    np.savetxt(f_rhom,np.array([0,0]).reshape(1,2),fmt='%s')                   
                    f_rhom.close()

                else:
                    
                    
                    unique_hosts_rand = np.unique(gdr).tolist()
                    subhalos = il.groupcat.loadSubhalos(basePath, j+1, fields = fSH)

                    ms_H = np.array(halos['GroupMassType'][:,4])*1e10/h
                    black_hole_mass_H = np.array(halos['GroupMassType'][:,5])*1e10/h
                    mstell_Model_H= np.array(ms_H[(np.unique(masked_desc)).astype(int)])
                    mstell_rModel_H = np.array(ms_H[(np.unique(gdr)).astype(int)])
                    mBH_H = ms_mBH(mstell_Model_H)
                    mBH_H = mBH_H[mBH_H>BHmass]
                    mBH_rH = ms_mBH(mstell_rModel_H)
                    mBH_rH = mBH_rH[mBH_rH>BHmass]
                    data_table = Table([np.repeat(j+1,len(mstell_Model_H),axis=None),np.unique(masked_desc),mstell_Model_H],names = ['#snap','#grnr','#mstell'],masked=True)
                    #save the stellar masses 
                    ascii.write(data_table,'stellar_mass_directory/stell_masses_snp'+str(j+1)+'_'+files[m]+'_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.csv', format='csv', overwrite=True)
                    data_tabler = Table([np.repeat(j+1,len(mstell_rModel_H),axis=None),np.unique(gdr),mstell_rModel_H],names = ['#snap','#grnr','#mstell'],masked=True)
                    ascii.write(data_tabler,'stellar_mass_directory/stell_massess_snp'+str(j+1)+'_'+files[m]+'_fpt'+seeding_str+'_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.csv', format='csv', overwrite=True)

                    #save nonzero mass densities
                    f_rhom = open('mass_density_directory/rhom_scalingrelation_mBH_ms_'+files[m]+'_fpt'+seeding_str+'_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.csv','a+')
                    np.savetxt(f_rhom,np.array([np.sum(mBH_H)/V,np.sum(mBH_rH)/V]).reshape(1,2),fmt='%s')
                    f_rhom.close()
                    
                    BHs_already_seeded.append(len(mBH_H))
                    BHs_already_seeded_rand.append(len(mBH_rH))
                        
                        
            if m in []: #indices for subhalo descendant table files only 
                
                get_desc= np.array(org_data_T[j]).astype(int)
                get_desc_prev = np.array(org_data_T[j-1]).astype(int)
                
                get_born_desc = get_desc[~np.isnan(get_desc)]
                get_born_desc_prev = get_desc_prev[~np.isnan(get_desc)]
                masked_desc = get_desc[(get_desc>=0)&(~np.isnan(get_desc))]
                print(len(masked_desc),len(get_born_desc))
                #masked_desc_prev = get_desc_prev[np.intersect1d(np.where(get_desc>=0)),(np.where()]
                random_SHs = np.random.randint(0, len(get_born_desc),np.round(int(seeding_percentage*len(get_born_desc)),0))
                random_SHs_list.extend(random_SHs)
                rs=np.array(np.unique(random_SHs_list), dtype = int)
                get_desc_rand = np.array(get_born_desc)[(rs)]
                get_desc_prev_rand = np.array(get_born_desc_prev)[(rs)]
                gdr = get_desc_rand[(get_desc_rand>=0)]
                pdr = get_desc_prev_rand[np.where(get_desc_rand>=0)]
                if len(masked_desc) == 0:
                    print("none")
                    nm.append(0)
                    BHs_already_seeded.append(0)
                    nm_rand.append(0)
                    BHs_already_seeded_rand.append(0)
                    f_rhom = open('mass_density_directory/rhom_scalingrelation_mBH_ms_'+files[m]+'_fpt'+seeding_str+'_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.csv','a+')
                    np.savetxt(f_rhom,np.array([0,0]).reshape(1,2),fmt='%s')                    
                    f_rhom.close()

                else:
                    
                    subhalos = il.groupcat.loadSubhalos(basePath, j+1, fields = fSH)
                    head = il.groupcat.loadHeader(basePath,j+1)
                    h = head['HubbleParam']
                    z = head['Redshift']

                    ms_SH = np.array(subhalos['SubhaloMassType'][:,4])*1e10/h
                    mstell_Model_SH= np.array(ms_SH[(np.unique(masked_desc)).astype(int)])
                    mstell_rModel_SH = np.array(ms_SH[(np.unique(gdr)).astype(int)])
                    mBH_SH = ms_mBH(mstell_Model_SH)
                    mBH_SH = mBH_SH[mBH_SH>BHmass]
                    mBH_rSH = ms_mBH(mstell_rModel_SH)
                    mBH_rSH = mBH_rSH[mBH_rSH>BHmass]
                    
                    BHs_already_seeded.append(len(mBH_SH))
                    BHs_already_seeded_rand.append(len(mBH_rSH))
                    f_rhom = open('mass_density_directory/rhom_scalingrelation_mBH_ms_'+files[m]+'_fpt'+seeding_str+'_'+res+'BHgr1e'+str(np.log10(BHmass))+'msun.csv','a+')
                    np.savetxt(f_rhom,np.array([np.sum(mBH_CSH)/V,np.sum(mBH_rCSH)/V]).reshape(1,2),fmt='%s')
                    f_rhom.close()
                    
         
        np.savetxt(fm,np.array(BHs_already_seeded),fmt='%s')
        np.savetxt(fmr,np.array(BHs_already_seeded_rand),fmt ='%s')

        fm.close()
        fmr.close()
        
