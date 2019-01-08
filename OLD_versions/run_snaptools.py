import eagle as E
import numpy as np
import snaptools as snaptools
from sys import exit
from sys import argv
import os

if 'overwrite' in argv:
    overwrite = True
else:
    overwrite = False

if 'addDM' in argv:
    addDM = True
else:
    addDM = False

if overwrite and addDM:
    print 'Please choose EITHER "overwrite" OR "addDM" OR no argument'
    exit()

if 'timing' in argv:
    timeflag = True
else:
    timeflag = False

##################################################################
# Available snapshots:
# ('redshift','snap')

snapdict = dict([('0.0','028_z000p000'),\
                 ('0.1','027_z000p101'),\
                 ('0.5','023_z000p503'),\
                 ('1.0','019_z001p004'),\
                 ('2.0','015_z002p012'),\
                 ('3.0','012_z003p017'),\
                 ('5.0','008_z005p037')])

##################################################################
# Choose what you want
simulation = 'L0100N1504'
run = 'REFERENCE'

#redshift_labels = ['0.5','1.0','3.0']

redshift_labels = ['3.0']

galtype = 'all' # 'dwarf', 'Lstar', 'group' or 'all'

particle_types = ['0','1','4']

save_path = '/data5/arijdav1/saved_regions/'

##################################################################

if addDM:
    particle_types = ['1','4']

for k in range(len(redshift_labels)):
    print 'Running halo isolator for ',simulation,'-',run
    print 'Redshift: ',redshift_labels[k]
    print 'Obtaining ',galtype,' galaxies...'

    snap = snapdict[redshift_labels[k]]
    
    if run == 'REFERENCE_ApogeeRun':
        sim = '/gal/simulations/EAGLE/L0025N0376/REFERENCE_ApogeeRun/data/'
    else:
        sim = '/data5/simulations/EAGLE/'+simulation+'/'+run+'/data/'

    # FOF quantities
    num_subs = np.array(E.readArray("SUBFIND_GROUP", sim, snap, "/FOF/NumOfSubhalos"))
    masslist = np.array(E.readArray("SUBFIND_GROUP",sim,snap,'FOF/Group_M_Crit200')[num_subs>0])*1e10
    first_subhalo = np.array(E.readArray("SUBFIND_GROUP",sim,snap,'FOF/FirstSubhaloID')[num_subs>0])
    r200s = np.array(E.readArray("SUBFIND_GROUP", sim, snap, "/FOF/Group_R_Crit200"))[num_subs>0]

    # SUBFIND quantities
    Mstar_30kpc = np.array(E.readArray("SUBFIND", sim, snap, "/Subhalo/ApertureMeasurements/Mass/030kpc"))[first_subhalo,4] * 1e10
    subfind_centres = np.array(E.readArray('SUBFIND',sim,snap,'Subhalo/CentreOfPotential'))[first_subhalo,:]
    subfind_subgroupnums = np.array(E.readArray('SUBFIND',sim,snap,'Subhalo/SubGroupNumber'))[first_subhalo]
    subfind_bulkvels = np.array(E.readArray('SUBFIND',sim,snap,'Subhalo/Velocity'))[first_subhalo,:]

    # Select your halos
    if galtype=='dwarf':
        mask = np.where((masslist>=10**11.3)&(masslist<=10**11.5))[0]
    elif galtype=='Lstar':
        mask = np.where((masslist>=10**11.9)&(masslist<=10**12.25))[0] # did go up to only 12.1 before
    elif galtype=='group':
        mask = np.where((masslist>=10**12.6)&(masslist<=10**13.2))[0]
    elif galtype=='all':
        mask = np.where((masslist>=10**11.5)|(Mstar_30kpc>=10**9.5))[0] # NEW CRITERION

    groupnums = mask+1
    r200s = r200s[mask]
    centres = subfind_centres[mask]
    subgroupnums = subfind_subgroupnums[mask]
    bulkvels= subfind_bulkvels[mask]

    print 'Scanning existing galaxy files...'
    nonexist_indices = []
    for c in range(len(centres)):
        savepath='/data5/arijdav1/saved_regions/'+simulation+'_'+run+'/'+snap+'/'+galtype+'/group'+str(groupnums[c])+'.hdf5'
        if not os.path.exists(savepath):
            nonexist_indices.append(c)

    if not overwrite and not addDM:
        if not nonexist_indices:
            print 'All selected galaxies have already been saved. To overwrite, enter "overwrite" as a command line argument'
            exit()

        print len(centres)-len(nonexist_indices),' of the ',len(centres),' groups selected have already been saved and will not be re-done.'
        print 'To overwrite, enter "overwrite" as a command line argument"'

        centres = centres[nonexist_indices]
        groupnums = groupnums[nonexist_indices]
        subgroupnums = subgroupnums[nonexist_indices]
        bulkvels = bulkvels[nonexist_indices]



    # Set the box size such that we can make profiles out to log(r/r200) = 0.5 of the biggest halo
    reg_sizes = 2 * r200s * 10**0.5

    print len(centres),' halos to do.'

    volume = snaptools.Snapshot_Volume(centres, reg_sizes, bulkvels, groupnums,subgroupnums,run=simulation,model=run,tag=snap,halotype=galtype,part_types=particle_types,timing_flag=timeflag)

    volume.save_individuals(path=save_path,overwrite=overwrite,addDM=addDM)


