import numpy as np
import h5py as h5
import os
from sys import exit
from sys import argv
from tqdm import tqdm
import eagle as E

def get_saved_halos(sim,run,tag,halotype):
    from glob import glob
    files = sorted(glob('/hpcdata7/arijdav1/saved_regions/'+sim+'_'+run+'/'+tag+'/'+halotype+'/group*.hdf5'))
    groupnums = []
    for f in range(len(files)):
        split1 = files[f].split('group',1)[1]
        groupnums.append(split1.split('.',1)[0])
    return np.sort(map(int,groupnums))

def ptype(parttype):
    ptdict = dict([('Gas','0'),('DM','1'),('Stars','4'),('BHs','5')])
    return ptdict[parttype]

def load_array(datafield):
    sim_location = '/hpcdata5/simulations/EAGLE/'+sim+'/'+run+'/data'
    return np.array(E.readArray(datatype,sim_location,tag,'/PartType'+ptype(parttype)+'/'+datafield))

################################################################
# Inputs
sim = 'L0025N0376'
run = 'REFERENCE'
tag='028_z000p000'
#tag = '023_z000p503'
datatype = 'SNAP'

halotype = 'all' # This field is for when I save my halos in mass bins, you can ignore it

parttype = 'Stars' # Particle type you want

datafield = 'StellarFormationTime' # Pick the data you want from the full snapshot

overwrite = False # If the data has already been downloaded to the file, overwrite it?
verbose = False # Print out useful things to the command line

################################################################

print 'Obtaining ',datafield

if datafield == 'SubGroupNumber':
    print 'Loading subgroup matching files...'

    if run == 'AGNdT8':
        sgn_file = h5.File('/hpcdata5/simulations/EAGLE/aux_data/'+sim+'/EagleVariation_AGNdT8p00/snap_grps/snapshot_grp_info_'+tag+'.hdf5','r')
    elif run == 'AGNdT9':
        sgn_file = h5.File('/hpcdata5/simulations/EAGLE/aux_data/'+sim+'/EagleVariation_AGNdT9p00/snap_grps/snapshot_grp_info_'+tag+'.hdf5','r')
    else:
        sgn_file = h5.File('/hpcdata5/simulations/EAGLE/aux_data/'+sim+'/'+run+'/snap_grps/snapshot_grp_info_'+tag+'.hdf5','r')
    eagle_pids = np.array(sgn_file['PartType0/ParticleIDs'])
    eagle_sgns = np.array(sgn_file['PartType0/SubGroupNumber'])

    match_sort = np.argsort(eagle_pids)
    sorted_eagle_data = eagle_sgns[match_sort]
    sorted_eagle_pids = eagle_pids[match_sort]

else:
    print 'Loading all EAGLE ParticleIDs...'
    eagle_pids = load_array('ParticleIDs').astype(int)
    sorted_eagle_pids = np.sort(eagle_pids)
    print 'Loading all EAGLE data for '+datafield
    sorted_eagle_data = load_array(datafield)[np.argsort(eagle_pids)]


saved_halos = get_saved_halos(sim,run,tag,halotype) # Get a list of the group numbers already downloaded. You'll need to alter the location of your saved files in the function at the top!

    
print 'Saving to snaptools halos...'
for n in tqdm(range(len(saved_halos))):

    try:
        # Specify the location of your downloaded halo (same as in get_saved_halos function above)
        f = h5.File('/hpcdata7/arijdav1/saved_regions/'+sim+'_'+run+'/'+tag+'/'+halotype+'/group'+str(saved_halos[n])+'.hdf5','r+')
    except IOError:
        print 'get_saved_halos has failed in this instance'
        exit()
    if verbose:
        print 'Loaded halo ',saved_halos[n]

    if not overwrite and datafield in f[parttype].keys():
        if verbose:
            print 'This data already exists in the halofile, overwrite mode is OFF.'
        continue

    if overwrite and datafield in f[parttype].keys():
        del f[parttype+'/'+datafield]

    halo_pids = np.array(f[parttype+'/ParticleIDs'],dtype=int)

    if verbose:
        print 'Getting '+datafield+' data for halo...'

    try:
        halo_newdata = sorted_eagle_data[np.searchsorted(sorted_eagle_pids,halo_pids)]
    except IndexError:
        continue

    if verbose:
        print 'Saving...'
    f[parttype].create_dataset(datafield,data=halo_newdata)

    



