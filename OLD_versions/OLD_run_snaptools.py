import eagle as E
import numpy as np
import snaptools as snaptools
from sys import exit

##################################################################

# Choose what you want
simulation = 'L0100N1504'
run = 'REFERENCE'
snap = '028_z000p000'



galtype = 'highmass' # 'dwarf', 'Lstar', 'group' or 'all' OR 'highmass' for original 100 box download




infotype = 'SUBFIND'
table = 'Subhalo'
# Pick by stellar mass or halo mass?
halo = True

##################################################################

# Simulation information
h = 0.6777
#z = 0.101
if simulation == 'L0100N1504':
	boxsize = 100
elif simulation == 'L0025N0752':
	boxsize = 25
elif simulation == 'L0025N0376':
	boxsize = 25
elif simulation == 'L0050N0752':
    boxsize = 50

##################################################################

vol = boxsize**3

sim = '/data5/simulations/EAGLE/'+simulation+'/'+run+'/data/'

print 'Simulation: ',simulation
print 'Model: ',run
print 'Halo type: ',galtype

if halo:             # Limits used in Oppenheimer+2016, but made narrower for stacking
    num_subs = np.array(E.readArray("SUBFIND_GROUP", sim, snap, "/FOF/NumOfSubhalos"))
    masslist = np.array(E.readArray("SUBFIND_GROUP",sim,snap,'FOF/Group_M_Crit200')[num_subs>0])*1e10

    
    if galtype=='dwarf':
        mask = np.where((masslist>=10**11.3)&(masslist<=10**11.5))[0]
    elif galtype=='Lstar':
        mask = np.where((masslist>=10**11.9)&(masslist<=10**12.25))[0] # this was to 12.1 before.
    elif galtype=='group':
        mask = np.where((masslist>=10**12.6)&(masslist<=10**13.2))[0]
    elif galtype=='all':
        mask = np.where((masslist>=10**11)&(masslist<=10**15))[0]
    elif galtype=='highmass':
        mask = np.where((masslist>=np.power(10.,11.5))&(masslist<=10**15))[0]
    gns = mask+1
    

    r200 = np.array(E.readArray("SUBFIND_GROUP", sim, snap, "/FOF/Group_R_Crit200"))[num_subs>0]
    r200 = r200[gns]

    subfind_centres = np.array(E.readArray(infotype,sim,snap,table+'/CentreOfPotential'))
    subfind_groupnums = np.array(E.readArray(infotype,sim,snap,table+'/GroupNumber'))
    subfind_subgroupnums = np.array(E.readArray(infotype,sim,snap,table+'/SubGroupNumber'))
    subfind_bulkvels = np.array(E.readArray(infotype,sim,snap,table+'/Velocity'))

    masktwo = []    
    for n, gn in enumerate(subfind_groupnums):
        if gn in gns and subfind_subgroupnums[n]==0:
            masktwo.append(n)

    centres = subfind_centres[masktwo]
    groupnums = subfind_groupnums[masktwo]
    subgroupnums = subfind_subgroupnums[masktwo]
    bulkvels= subfind_bulkvels[masktwo]


else:
    massarr = E.readArray(infotype,sim,snap,table+'/ApertureMeasurements/Mass/030kpc')*1e10
    massarr = np.array(massarr)
    masslist = massarr[:,4]

    centres = E.readArray(infotype,sim,snap,table+'/CentreOfPotential')
    groupnums = E.readArray(infotype,sim,snap,table+'/GroupNumber')
    subgroupnums = E.readArray(infotype,sim,snap,table+'/SubGroupNumber')
    bulkvels = E.readArray(infotype,sim,snap,table+'/Velocity')

    maskone = np.where(subgroupnums==0)[0] # ONLY CENTRALS - IMPORTANT FOR FINDING R200

    masslist = masslist[maskone]
    centres = centres[maskone,:]
    groupnums = groupnums[maskone]
    subgroupnums = subgroupnums[maskone]
    bulkvels = bulkvels[maskone,:]

    masktwo = np.where((masslist>=1e10)&(masslist<=2e10))[0]

    masslist = masslist[masktwo]
    centres = centres[masktwo,:]
    groupnums = groupnums[masktwo]
    subgroupnums = subgroupnums[masktwo]
    bulkvels = bulkvels[masktwo,:]

    #print 'Found ',len(masslist),' galaxies in this mass range'


# Set the box size such that we can make profiles out to log(r/r200) = 0.5 of the biggest halo
reg_sizes = 2 * r200 * 10**0.5

print len(centres)

volume = snaptools.Snapshot_Volume(centres, reg_sizes, bulkvels, groupnums,subgroupnums,run=simulation,model=run,tag=snap,halotype=galtype)
volume.save_individuals()


