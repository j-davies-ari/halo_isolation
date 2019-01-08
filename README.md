# halo_isolation

halo_isolation is a module for effectively "chopping out galaxies" from EAGLE snapshots. Given a set of FOF group numbers, centres of potential and bulk halo velocities, it will create new HDF5 files containing only the particles within a distance of your choosing from each halo's centre. Pre-processing EAGLE snapshots in this way can greatly speed up and/or simplify your code if you are working with individual galaxies, as you won't need to load in and manipulate entire EAGLE snapshots every time you want to run your code. The new HDF5 files have a very similar structure to EAGLE snapshots and contain most of the important particle quantities - you can always add more quantities within the module's code.

The file "run_halo_isolator.py" is what I use to run the module - take a look at it and the included comments to see how to use the module.

There are only a few dependencies:

- The standard "readEagle" module used for reading EAGLE snapshots - this can be found on the EAGLE wiki under "Reading snapshots and group files in python"

- The rather confusingly named "read_eagle" module, which is not the same module as above - it uses the hash tables to very quickly load in small regions of the snapshots and does most of the heavy lifting in this module. It can also be obtained from the EAGLE wiki under "Read routines for PH key sorted snapshots".

- "h5py", which is necessary for saving the hdf5 output

- "tqdm", which is a lovely little module for adding a nice, helpful progress bar to your terminal during for loops. 

- "metadata", which contains functions for finding EAGLE snapshot tags at a given redshift, as well as several helpful constants and units for the simulations. I have included this in the repository.
