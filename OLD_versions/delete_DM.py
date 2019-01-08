import numpy as np
import h5py as h5
from sys import exit
import eagle as E
from tqdm import tqdm

sim = 'L0100N1504'
run = 'REFERENCE'
snap = '028_z000p000'
halotype = 'all'

# how many galaxies are there to do?
simloc = '/data5/simulations/EAGLE/'+sim+'/'+run+'/data/'
num_subs = np.array(E.readArray("SUBFIND_GROUP", simloc, snap, "/FOF/NumOfSubhalos"))
masslist = np.array(E.readArray("SUBFIND_GROUP",simloc,snap,'FOF/Group_M_Crit200')[num_subs>0])*1e10
first_subhalo = np.array(E.readArray("SUBFIND_GROUP",simloc,snap,'FOF/FirstSubhaloID')[num_subs>0])
Mstar_30kpc = np.array(E.readArray("SUBFIND", simloc, snap, "/Subhalo/ApertureMeasurements/Mass/030kpc"))[first_subhalo,4] * 1e10
gns = np.arange(0,len(masslist)) + 1
gns = gns[(masslist>np.power(10.,11.5))|(Mstar_30kpc>np.power(10.,9.5))]

print 'Are you sure you want to delete the DM data for ',sim,'-',run,'-',snap,'? (y/n)'
answer = raw_input('---->')
if answer != 'y':
    exit()

for k in tqdm(range(len(gns))):
    n = gns[k]

    try:
        f = h5.File('/data5/arijdav1/saved_regions/'+sim+'_'+run+'/028_z000p000/'+halotype+'/group'+str(n)+'.hdf5','r+')
    except IOError:
        continue

    del f['DarkMatter/Coordinates']
    del f['DarkMatter/Velocity']
    del f['DarkMatter/Mass']
    del f['DarkMatter/PID']
    del f['DarkMatter/GroupNumber']

    f.close()
