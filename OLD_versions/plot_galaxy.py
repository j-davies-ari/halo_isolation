import numpy as np
import h5py as h5
from sys import exit
from plot_tools import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

for n in range(100):
	print 'Loading data'
	try:
		data = h5.File('/home/arijdav1/Desktop/L0025N0752_REFERENCE/028_z000p000/group'+str(n)+'.hdf5','r')
	except IOError:
		continue
	print 'Loaded.'

	star_coords = np.array(data['Stars/Coordinates'])
	gas_coords = np.array(data['Gas/Coordinates'])
	print len(star_coords[:,0]),' star particles to plot'
	print len(gas_coords[:,0]),' gas particles to plot'

	# Convert to kpc
	star_coords *= 1e3
	gas_coords *= 1e3

	xedges = np.linspace(-10,10,num=200)
	yedges = np.linspace(-10,10,num=200)

	starmap,xedges,yedges = np.histogram2d(star_coords[:,1],star_coords[:,0],bins=(xedges,yedges))
	gasmap,xedges,yedges = np.histogram2d(gas_coords[:,1],gas_coords[:,0],bins=(xedges,yedges))

	xcents = get_bincentres(xedges)
	ycents = get_bincentres(yedges)
	X, Y = np.meshgrid(xedges, yedges)


	'''
	fig = plt.figure()
	ax = Axes3D(fig)
	ax.scatter(star_coords[:,0],star_coords[:,1],star_coords[:,2],c='k')
	plt.show()
	'''

	fig = plt.figure(figsize=(10,10))
	'''
	ax = fig.add_subplot(111)
	ax.pcolormesh(X, Y, starmap)
	ax.set_aspect('equal')
	ax.set_xlim(-10,10)
	ax.set_ylim(-10,10)
	'''
	#im2 = plt.imshow(starmap,interpolation='bilinear',cmap='bone')
	#plt.hold(True)
	im1 = plt.imshow(gasmap,interpolation='bilinear',cmap='hot')
	plt.scatter(star_coords[:,1],star_coords[:,0],c='b')
	plt.show()
