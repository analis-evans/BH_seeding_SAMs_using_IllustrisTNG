import sys
import numpy as np
import h5py
import arepo_package
import scipy.interpolate
import illustris_python as il
import illustris_python.sublink
import time
from mpi4py import MPI
import os
import psutil
import multiprocessing

comm = MPI.COMM_WORLD
process = psutil.Process(os.getpid())
start = time.process_time()
rank = MPI.COMM_WORLD.rank
nproc = MPI.COMM_WORLD.Get_size()   #Size of communicator                       
dp = nproc
iproc = MPI.COMM_WORLD.Get_rank()   #Ranks in communicator

n1 = 
n2 = 
n3 = 
n4 = 
#Modify file list with metallicity cut strings ('H'denotes halo, 'SH' denotes 'subhalo','cSH' denotes 'central subhalo', 'zpt1' denotes Z < 0.1, 'sfr' denotes star-forming)
headers=['mg'+str(n1)+'e'+str(n3)+'_mH'+str(n2)+'e'+str(n4)+'_zpt1',
         'mg'+str(n1)+'e'+str(n3)+'_mH'+str(n2)+'e'+str(n4)+'_zpt001',
         'mg'+str(n1)+'e'+str(n3)+'_mH'+str(n2)+'e'+str(n4)+'_zpt1sfr',
         'mg'+str(n1)+'e'+str(n3)+'_mH'+str(n2)+'e'+str(n4)+'_zpt001sfr',
         'mg1e7'+str(n1)+'e'+str(n3)+'_mSH1e8'+str(n2)+'e'+str(n4)+'_zpt1',
         'mg'+str(n1)+'e'+str(n3)+'_mSH'+str(n2)+'e'+str(n4)+'_zpt001',
         'mg'+str(n1)+'e'+str(n3)+'_mSH'+str(n2)+'e'+str(n4)+'_zpt1sfr',
         'mg'+str(n1)+'e'+str(n3)+'_mSH'+str(n2)+'e'+str(n4)+'_zpt001sfr',
         'mg'+str(n1)+'e'+str(n3)+'_mcSH'+str(n2)+'e'+str(n4)+'_zpt1',
         'mg'+str(n1)+'e'+str(n3)+'_mcSH'+str(n2)+'e'+str(n4)+'_zpt001',
         'mg'+str(n1)+'e'+str(n3)+'_mcSH'+str(n2)+'e'+str(n4)+'_zpt1sfr',
         'mg'+str(n1)+'e'+str(n3)+'_mcSH'+str(n2)+'e'+str(n4)+'_zpt001sfr',
         'TNGcut_mH5e10msolperh'
]

path_to_simulation='path_to_simulation'

run= 'L35n2160TNG'                           

basePath=path_to_simulation + run + '/output/'

path_to_TNG_trees=path_to_simulation + run  +'/postprocessing/trees/SubLink/'

res =
m = 
n = int(sys.argv[1]) #corresponds to argument $1 from the mpiexec command in slurm file; the number of cores must be set to the number of subfiles in the TNG tree run folder; the loop in the *.sbatch file must then be set to (0 1 seq-1) 

snapshots_at_unique_IDs = np.array(np.loadtxt('path/'+headers[m]+res+'subhalos.txt',usecols = 0),order = 'K')                                                                             
subhalo_indices = np.array(np.loadtxt('path/'+headers[m]+res+ 'subhalos.txt', usecols = 1),order = 'K')

seq = 8 #set this to the number of parts you want to separate the subhalo and snapshot lists into for parallel processing (must be a number reasonably smaller than the length of the list) ;  the loop in the *.sbatch file must then be set to (0 1 seq-1)
list_by_dp = int(np.round(len(subhalo_indices)/(seq),0))
file_desc = open('desc_'+headers[m]+res+'_sysargv'+str(n)+'.csv','a+')
if n==(seq-1):
    subhalo_indices = subhalo_indices[list_by_dp*(n):] 
    snapshots_at_unique_IDs = snapshots_at_unique_IDs[list_by_dp*(n):]

if n<=(seq-2):
    subhalo_indices = subhalo_indices[list_by_dp*(n):list_by_dp*(n+1)]
    snapshots_at_unique_IDs = snapshots_at_unique_IDs[list_by_dp*(n):list_by_dp*(n+1)]

SubfindID, SnapNum, SubhaloID, DescendantID, TreeID, new_SH_list,new_snap_list=arepo_package.load_tree_for_descendants(basePath,TNG=1,FileNum=iproc, snap_list=snapshots_at_unique_IDs,subhalo_list=subhalo_indices,path_to_TNG_trees=path_to_TNG_trees)

for k in range(0,len(new_SH_list)):
    descendant_indices1, snapshot_indices1 = arepo_package.get_sublink_descendants(subhalo_index=new_SH_list[k],output_snapshot= new_snap_list[k])

    if (snapshot_indices1[-1]!= 99):
        end_list = np.repeat(-1000, 99-snapshot_indices1[-1])
        descendant_indices1.extend(end_list)
                                                                            
    xind1 = np.in1d(np.arange(snapshot_indices1[0],snapshot_indices1[-1]+1,1), snapshot_indices1)
    x_ind1 = np.where(xind1*1 == 0)[0]
    if len(x_ind1) > 0:
        for q in x_ind1.tolist():
            descendant_indices1.insert(q,-1000) #place holder for descendants with missing/unidentified links in the TNG merger tree 
    nan_array1 = np.repeat(np.nan, 99-len(descendant_indices1))
    data1 = np.concatenate((nan_array1, descendant_indices1), axis=None)
    print(data1,flush=True)
    file_desc = np.savetxt(file_desc, data1.reshape(1,99), fmt ='%s')

MPI.Finalize()


