import numpy as np
import h5py as h5
from plot_tools import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from matplotlib import cm

data = h5.File('/data5/simulations/EAGLE/L0025N0376/REFERENCE/data/snapshot_028_z000p000/snap_028_z000p000.0.hdf5','r')


coords = np.array(data['PartType4/Coordinates'])
groupno = np.array(data['PartType4/GroupNumber'])
temp = np.array(data['PartType4/MaximumTemperature'])
subgroupno = np.array(data['PartType4/SubGroupNumber'])
holes = np.array(data['PartType5/Coordinates'])
hole_groupno = np.array(data['PartType5/GroupNumber'])


#######################################################################################
# Find out how many particles are associated with each group number

sortgroup = np.sort(groupno)
unique_gnos = sortgroup[unique_entries(sortgroup)]
out = np.zeros((len(unique_gnos),2))
out[:,0]=unique_gnos
for i in range(len(unique_gnos)):
	wh = np.where(groupno==unique_gnos[i])
	out[i,1] = len(wh[0])
print out.astype(int)

#######################################################################################




group = 12

testgal_coords = coords[groupno==group]
testgal_temp = temp[groupno==group]
testgal_subgroupno = subgroupno[groupno==group]
testgal_holes = holes[hole_groupno==group]

print len(testgal_coords),' particles in this group'
'''
print testgal_coords
print testgal_subgroupno
'''

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(testgal_coords[:,0],testgal_coords[:,1],testgal_coords[:,2],c=testgal_temp,cmap = 'Reds')
#ax.scatter(testgal_holes[:,0],testgal_holes[:,1],testgal_holes[:,2],c='cyan',marker='o',s=10000)
plt.show()










