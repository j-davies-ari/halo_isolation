import numpy as np
import h5py as h5
import os
from sys import exit
from sys import argv
from tqdm import tqdm
import eagle as E

def get_saved_halos(sim,run,tag,halotype):
    from glob import glob
    files = sorted(glob('/data5/arijdav1/saved_regions/'+sim+'_'+run+'/'+tag+'/'+halotype+'/group*.hdf5'))
    groupnums = []
    for f in range(len(files)):
        split1 = files[f].split('group',1)[1]
        groupnums.append(split1.split('.',1)[0])
    return np.sort(map(int,groupnums))

# Inputs
simulation = 'L0100N1504'
run = 'REFERENCE'
tag='028_z000p000'

halotype = 'all'

sim = '/data5/simulations/EAGLE/'+simulation+'/'+run+'/data/'

saved_halos = get_saved_halos(simulation,run,tag,halotype) # Get a list of the group numbers already downloaded. You'll need to alter the location of your saved files in the function at the top!

num_subs = np.array(E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/NumOfSubhalos"))
r200s = np.array(E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/Group_R_Crit200"))[num_subs>0]

r200s = r200s[saved_halos-1]

print 'Saving the new r200s'
for n in tqdm(range(len(saved_halos))):
    try:
        # Specify the location of your downloaded halo (same as in get_saved_halos function above)
        f = h5.File('/data5/arijdav1/saved_regions/'+simulation+'_'+run+'/028_z000p000/'+halotype+'/group'+str(saved_halos[n])+'.hdf5','r+')
    except IOError:
        print 'get_saved_halos has failed in this instance'
        exit()

    del f['Volume/r200']
    f['Volume'].create_dataset('r200',data=r200s[n])

    f.close()
