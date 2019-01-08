import eagle as E
import numpy as np
from sys import exit
from sys import argv
import os
import halo_isolator
from tqdm import tqdm
import metadata

if 'overwrite' in argv:
    overwrite = True
else:
    overwrite = False

##################################################################
# Choose which simulation you want
simulation = 'L0050N0752'
run = 'EagleVariation_NOAGN'

# A list of redshifts that you want to save halos from
redshifts = [0.,]

# this label is used when saving the halos - you can remove this functionality, see below.
halotype = 'all' # 'dwarf', 'Lstar', 'group' or 'all'

# Change these directories to suit your system
save_path = '/hpcdata7/arijdav1/saved_regions/' # the location for the output hdf5 files

data_location = '/hpcdata5/simulations/EAGLE/' # point to the EAGLE data - this is the LJMU directory

##################################################################

# Loops over redshifts
for k in range(len(redshifts)):
    print 'Running halo isolator for ',simulation,'-',run
    print 'Redshift: ',redshifts[k]
    print 'Obtaining ',halotype,' galaxies...'

    # This uses my little 'metadata' module to find the snapshot tag for your chosen redshift
    tag =  metadata.snapnum_search(redshifts[k],returntag=True)
    sim = data_location+simulation+'/'+run+'/data/'

    #tag = '999_z000p000'

    ###############################################################################################
    '''
    Put some code here to generate lists of FOF group numbers, their centres of potential and
    their bulk halo velocities, which the module subtracts off. I select my halos by mass,
    but you could always just input a list of a few halos of interest. The only quirk is that
    the save() function takes an input called 'halotype' which I use to separate my mass bins,
    you can just remove this functionality from save() in the module if you like - uncomment the
    'vol_dir' line which makes the directory structure more EAGLE-like and remove 'halotype' from
    the function arguments.
    '''

    # FOF quantities
    num_subs = np.array(E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/NumOfSubhalos"))
    masslist = np.array(E.readArray("SUBFIND_GROUP",sim,tag,'FOF/Group_M_Crit200'))*1e10
    first_subhalo = np.array(E.readArray("SUBFIND_GROUP",sim,tag,'FOF/FirstSubhaloID'))
    r200s = np.array(E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/Group_R_Crit200"))
    groupnums = np.arange(len(num_subs)) + 1

    # SUBFIND quantities

    Mstar_30kpc = np.array(E.readArray("SUBFIND", sim, tag, "/Subhalo/ApertureMeasurements/Mass/030kpc"))[:,4] * 1e10
    max_subhalo = len(Mstar_30kpc)
    # If the final FOF group is empty, it can be assigned a FirstSubhaloID which is out of range
    first_subhalo[first_subhalo==max_subhalo] -= 1 # this line fixes the issue
    Mstar_30kpc = Mstar_30kpc[first_subhalo]
    subfind_centres = np.array(E.readArray('SUBFIND',sim,tag,'Subhalo/CentreOfPotential'))[first_subhalo,:]
    subfind_bulkvels = np.array(E.readArray('SUBFIND',sim,tag,'Subhalo/Velocity'))[first_subhalo,:]

    # Remove empty FOF groups
    nonempty_groups = np.where(num_subs>0)[0]
    masslist = masslist[nonempty_groups]
    Mstar_30kpc = Mstar_30kpc[nonempty_groups]
    r200s = r200s[nonempty_groups]
    subfind_centres = subfind_centres[nonempty_groups]
    subfind_bulkvels = subfind_bulkvels[nonempty_groups]
    groupnums = groupnums[nonempty_groups]

    # Select your halos by mass
    if halotype=='dwarf':
        mask = np.where((masslist>=10**11.3)&(masslist<=10**11.5))[0]
    elif halotype=='Lstar':
        mask = np.where((masslist>=10**11.9)&(masslist<=10**12.25))[0] # did go up to only 12.1 before
    elif halotype=='group':
        mask = np.where((masslist>=10**12.6)&(masslist<=10**13.2))[0]
    elif halotype == 'lowmass':
        mask = np.where((masslist>=10**11.)&(masslist<=10**11.5))[0]
    elif halotype=='all':
        mask = np.where((masslist>=10**11.5)|(Mstar_30kpc>=10**9.5))[0] # NEW CRITERION

    # These are the three arrays you NEED
    groupnums = groupnums[mask]
    centres = subfind_centres[mask]
    bulkvels= subfind_bulkvels[mask]

    # I use r200 of each halo to set my region size to cut out
    r200s = r200s[mask]

    ###############################################################################################

    '''
    This bit will check your save directory to see if any of your selected galaxies have already
    been saved - if you haven't entered 'overwrite' as a command line argument it will skip over
    any halos that already exist.
    '''

    print 'Scanning existing galaxy files...'
    nonexist_indices = []
    for c in range(len(centres)):
        savepath=save_path+simulation+'_'+run+'/'+tag+'/'+halotype+'/group'+str(groupnums[c])+'.hdf5'
        if not os.path.exists(savepath):
            nonexist_indices.append(c)

    if not overwrite:
        if not nonexist_indices:
            print 'All selected galaxies have already been saved. To overwrite, enter "overwrite" as a command line argument'
            exit()

        print len(centres)-len(nonexist_indices),' of the ',len(centres),' groups selected have already been saved and will not be re-done.'
        print 'To overwrite, enter "overwrite" as a command line argument"'

        centres = centres[nonexist_indices]
        groupnums = groupnums[nonexist_indices]
        bulkvels = bulkvels[nonexist_indices]

    ###############################################################################################

    # Here's where we do the work and run the module.

    print len(centres),' halos to do.'

    # Create a list of region sizes to chop out here which is the same length as the number of halos you want
    # I go out to a multiple of r200 for each galaxy
    reg_sizes = 2 * r200s * 10**0.5

    # Initialise the module for the simulation output we want
    box = halo_isolator.snapshot(sim=simulation,run=run,tag=tag,data_location=data_location)

    # Loop over the FOF halos. The 'tqdm' module gives you a nice progress bar!
    for i in tqdm(range(len(centres))):
        # Most of the work is done here!
        chunk = box.read_chunk(groupnums[i],centres[i],bulkvels[i],reg_sizes[i])
        # Generates a new .hdf5 file in your chosen location for each halo.
        box.save(chunk,halotype,path=save_path,overwrite=overwrite)



