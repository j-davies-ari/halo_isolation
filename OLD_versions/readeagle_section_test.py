import read_eagle
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

# Name of one file from the snapshot
snapfile = "/data5/simulations/EAGLE/L0100N1504/REFERENCE/data/snapshot_028_z000p000/snap_028_z000p000.0.hdf5"

print 'Opening snapshot'
# Open snapshot
snap = read_eagle.EagleSnapshot(snapfile)

# Select region of interest
snap.select_region(0, 10, 0, 10, 0, 10)

print 'Splitting process'
# Split selection between processors
# This assigns an equal number of hash cells to each processor.
snap.split_selection(comm_rank, comm_size)

print 'Reading data'
# Read data - each processor will receive a part of the selected region
density = snap.read_dataset(0, "Density")
temperature = snap.read_dataset(0, "Temperature")
pos = snap.read_dataset(0, "Coordinates")
ids = snap.read_dataset(0, "ParticleIDs")

print 'Making plot'
plt.figure()
plt.scatter(np.log10(density),np.log10(temperature))

plt.show()