import numpy as np
import os
import sys
import eagle as E
import read_eagle as read
import h5py
import copy



class snapshot(object):

    def __init__(self,sim = 'L0100N1504',
                        run = 'REFERENCE',
                        tag='028_z000p000',
                        pdata_type = 'SNAPSHOT',
                        data_location = '/hpcdata5/simulations/EAGLE/'):

        # Initialises everything and loads in information about the entire simulation volume

        if run == 'EagleVariation_NOAGN':
            self.incBH = False
        else:
            self.incBH = True

        self.sim = sim
        self.run = run
        self.tag = tag

        sim_path = data_location + sim + '/' + run + '/data/'

        # Name of one file from the snapshot
        self.snapfile = sim_path+'snapshot_'+tag+'/snap_'+tag+'.0.hdf5'

        # Get volume information
        boxsize = E.readAttribute(pdata_type, sim_path, tag, "/Header/BoxSize")
        self.h = E.readAttribute(pdata_type, sim_path, tag, "/Header/HubbleParam")
        self.Omega0 = E.readAttribute(pdata_type, sim_path, tag, "/Header/Omega0")
        self.OmegaLambda = E.readAttribute(pdata_type, sim_path, tag, "/Header/OmegaLambda")
        self.OmegaBaryon = E.readAttribute(pdata_type, sim_path, tag, "/Header/OmegaBaryon")
        self.a_0 = E.readAttribute(pdata_type, sim_path, tag, "/Header/ExpansionFactor")
        self.physical_boxsize = boxsize * self.a_0/self.h
        self.masstable = E.readAttribute(pdata_type, sim_path, tag, "/Header/MassTable") / self.h
        self.r200s = E.readArray("SUBFIND_GROUP", sim_path, tag, "/FOF/Group_R_Crit200")



    def read_chunk(self,groupnumber,centre,bulkvel,regionsize):

        '''
        This function uses the read_eagle module to very quickly load in a chunk of a snapshot
        around a particular FOF group with side length equal to the region size you want to save.
        It centres the co-ordinates on the centre of potential of the group (which you have input),
        and wraps the periodic box. It then masks this chunk to a spherical region with a diameter
        equal to the region size. Using this mask, the module loads in many other properties for the
        particles in the region and generates a dictionary to store them. This is done for each
        particle type, then all the dictionaries are wrapped up in another one and output at the end.
        '''

        #print 'Reading gas particle data'
        ptype = 0
        gas_dict = {}

        code_centre = centre * self.h/self.a_0 # convert to h-less comoving code units
        code_regionsize = regionsize * self.h/self.a_0
        code_boxsize = self.physical_boxsize * self.h/self.a_0

        r200 = self.r200s[groupnumber-1]

        # Open snapshot
        snap = read.EagleSnapshot(self.snapfile)

        # Select region of interest
        snap.select_region(code_centre[0] - code_regionsize / 2.,
                            code_centre[0] + code_regionsize / 2.,
                            code_centre[1] - code_regionsize / 2.,
                            code_centre[1] + code_regionsize / 2.,
                            code_centre[2] - code_regionsize / 2.,
                            code_centre[2] + code_regionsize / 2.,)


        pos = snap.read_dataset(0,'Coordinates') # in code units

        pos -= code_centre # centre on galaxy

        # Wrap box
        pos[pos[:,0]<(-1.*code_boxsize/2.),0] += code_boxsize
        pos[pos[:,1]<(-1.*code_boxsize/2.),1] += code_boxsize
        pos[pos[:,2]<(-1.*code_boxsize/2.),2] += code_boxsize
        pos[pos[:,0]>code_boxsize/2.,0] -= code_boxsize
        pos[pos[:,1]>code_boxsize/2.,1] -= code_boxsize
        pos[pos[:,2]>code_boxsize/2.,2] -= code_boxsize

        # Convert to physical units
        pos *= self.a_0/self.h

        # Mask to region size
        rmax = regionsize/2.
        r2 = np.einsum('...j,...j->...',pos,pos) # get the radii from the centre
        mask = np.where(r2<rmax**2)[0] # make the mask        

        gas_dict['Coordinates'] = pos[mask]

        # Load in everything else. If you want any other quantities, simply add a similar line here
        # Don't forget to convert with factors of a and h as done here!
        gas_dict['Velocity'] = snap.read_dataset(ptype, "Velocity")[mask,:]*np.sqrt(self.a_0)  - bulkvel # subtract off halo velocity
        gas_dict['Mass'] = snap.read_dataset(ptype, "Mass")[mask] /self.h
        gas_dict['Density'] = snap.read_dataset(ptype, "Density")[mask] * self.h**2/self.a_0**3
        gas_dict['Temperature'] = snap.read_dataset(ptype, "Temperature")[mask]
        gas_dict['ParticleIDs'] = snap.read_dataset(ptype, "ParticleIDs")[mask]
        gas_dict['GroupNumber'] = snap.read_dataset(ptype, "GroupNumber")[mask]
        gas_dict['StarFormationRate'] = snap.read_dataset(ptype, "StarFormationRate")[mask]
        gas_dict['Metallicity'] = snap.read_dataset(ptype, "Metallicity")[mask]
        gas_dict['SmoothedMetallicity'] = snap.read_dataset(ptype, "SmoothedMetallicity")[mask]
        gas_dict['SmoothingLength'] = snap.read_dataset(ptype, "SmoothingLength")[mask] * self.a_0/self.h
        gas_dict['OnEquationOfState'] = snap.read_dataset(ptype, "OnEquationOfState")[mask]
        gas_dict['MaximumTemperature'] = snap.read_dataset(ptype, "MaximumTemperature")[mask]
        gas_dict['AExpMaximumTemperature'] = snap.read_dataset(ptype, "AExpMaximumTemperature")[mask]
        gas_dict['InternalEnergy'] = snap.read_dataset(ptype, "InternalEnergy")[mask]

        # I wrap up all the smoothed abundances into one big array, you can save them individually if you want.
        # Just create individaul entries in gas_dict for each element instead
        H = snap.read_dataset(ptype, "SmoothedElementAbundance/Hydrogen")[mask]
        He = snap.read_dataset(ptype, "SmoothedElementAbundance/Helium")[mask]
        C = snap.read_dataset(ptype, "SmoothedElementAbundance/Carbon")[mask]
        N = snap.read_dataset(ptype, "SmoothedElementAbundance/Nitrogen")[mask]
        O = snap.read_dataset(ptype, "SmoothedElementAbundance/Oxygen")[mask]
        Ne = snap.read_dataset(ptype, "SmoothedElementAbundance/Neon")[mask]
        Mg = snap.read_dataset(ptype, "SmoothedElementAbundance/Magnesium")[mask]
        Si = snap.read_dataset(ptype, "SmoothedElementAbundance/Silicon")[mask]
        S = Si*0.6054160
        Ca = Si*0.0941736
        Fe = snap.read_dataset(ptype, "SmoothedElementAbundance/Iron")[mask]
        gas_dict['Abundances'] = np.dstack((H,He,C,N,O,Ne,Mg,Si,S,Ca,Fe))[0]





        #print 'Reading dark matter particle data'
        ptype = 1
        DM_dict = {}

        pos = snap.read_dataset(1,'Coordinates') # in code units
        pos -= code_centre # centre on galaxy

        # Wrap box
        pos[pos[:,0]<(-1.*code_boxsize/2.),0] += code_boxsize
        pos[pos[:,1]<(-1.*code_boxsize/2.),1] += code_boxsize
        pos[pos[:,2]<(-1.*code_boxsize/2.),2] += code_boxsize
        pos[pos[:,0]>code_boxsize/2.,0] -= code_boxsize
        pos[pos[:,1]>code_boxsize/2.,1] -= code_boxsize
        pos[pos[:,2]>code_boxsize/2.,2] -= code_boxsize

        # Convert to physical units
        pos *= self.a_0/self.h

        # Mask to rmax
        r2 = np.einsum('...j,...j->...',pos,pos) # get the radii from the centre
        mask = np.where(r2<rmax**2)[0] # make the mask        

        DM_dict['Coordinates'] = pos[mask]

        # Load in everything else
        DM_dict['Velocity'] = snap.read_dataset(ptype, "Velocity")[mask,:]*np.sqrt(self.a_0)  - bulkvel
        DM_dict['Mass'] = np.ones(len(mask))*self.masstable[1]
        DM_dict['ParticleIDs'] = snap.read_dataset(ptype, "ParticleIDs")[mask]
        DM_dict['GroupNumber'] = snap.read_dataset(ptype, "GroupNumber")[mask]



        #print 'Reading star particle data'
        ptype = 4
        star_dict = {}

        pos = snap.read_dataset(4,'Coordinates') # in code units
        pos -= code_centre # centre on galaxy

        # Wrap box
        pos[pos[:,0]<(-1.*code_boxsize/2.),0] += code_boxsize
        pos[pos[:,1]<(-1.*code_boxsize/2.),1] += code_boxsize
        pos[pos[:,2]<(-1.*code_boxsize/2.),2] += code_boxsize
        pos[pos[:,0]>code_boxsize/2.,0] -= code_boxsize
        pos[pos[:,1]>code_boxsize/2.,1] -= code_boxsize
        pos[pos[:,2]>code_boxsize/2.,2] -= code_boxsize

        # Convert to physical units
        pos *= self.a_0/self.h

        # Mask to rmax
        r2 = np.einsum('...j,...j->...',pos,pos) # get the radii from the centre
        mask = np.where(r2<rmax**2)[0] # make the mask        

        star_dict['Coordinates'] = pos[mask]


        # Load in everything else
        star_dict['Velocity'] = snap.read_dataset(ptype, "Velocity")[mask,:]*np.sqrt(self.a_0)  - bulkvel
        star_dict['Mass'] = snap.read_dataset(ptype, "Mass")[mask] /self.h
        star_dict['Density'] = snap.read_dataset(ptype, "BirthDensity")[mask] * self.h**2/self.a_0**3
        star_dict['ParticleIDs'] = snap.read_dataset(ptype, "ParticleIDs")[mask]
        star_dict['GroupNumber'] = snap.read_dataset(ptype, "GroupNumber")[mask]
        star_dict['Metallicity'] = snap.read_dataset(ptype, "Metallicity")[mask]
        star_dict['SmoothingLength'] = snap.read_dataset(ptype, "SmoothingLength")[mask] * self.a_0/self.h

        H = snap.read_dataset(ptype, "SmoothedElementAbundance/Hydrogen")[mask]
        He = snap.read_dataset(ptype, "SmoothedElementAbundance/Helium")[mask]
        C = snap.read_dataset(ptype, "SmoothedElementAbundance/Carbon")[mask]
        N = snap.read_dataset(ptype, "SmoothedElementAbundance/Nitrogen")[mask]
        O = snap.read_dataset(ptype, "SmoothedElementAbundance/Oxygen")[mask]
        Ne = snap.read_dataset(ptype, "SmoothedElementAbundance/Neon")[mask]
        Mg = snap.read_dataset(ptype, "SmoothedElementAbundance/Magnesium")[mask]
        Si = snap.read_dataset(ptype, "SmoothedElementAbundance/Silicon")[mask]
        S = Si*0.6054160
        Ca = Si*0.0941736
        Fe = snap.read_dataset(ptype, "SmoothedElementAbundance/Iron")[mask]
        star_dict['Abundances'] = np.dstack((H,He,C,N,O,Ne,Mg,Si,S,Ca,Fe))[0]




        if self.incBH:
            #print 'Reading black hole particle data'
            ptype = 5
            BH_dict = {}

            pos = snap.read_dataset(5,'Coordinates') # in code units

            pos -= code_centre # centre on galaxy

            # Wrap box
            pos[pos[:,0]<(-1.*code_boxsize/2.),0] += code_boxsize
            pos[pos[:,1]<(-1.*code_boxsize/2.),1] += code_boxsize
            pos[pos[:,2]<(-1.*code_boxsize/2.),2] += code_boxsize
            pos[pos[:,0]>code_boxsize/2.,0] -= code_boxsize
            pos[pos[:,1]>code_boxsize/2.,1] -= code_boxsize
            pos[pos[:,2]>code_boxsize/2.,2] -= code_boxsize

            # Convert to physical units
            pos *= self.a_0/self.h

            # Mask to rmax
            r2 = np.einsum('...j,...j->...',pos,pos) # get the radii from the centre
            mask = np.where(r2<rmax**2)[0] # make the mask        

            BH_dict['Coordinates'] = pos[mask]

            # Load in everything else
            BH_dict['Velocity'] = snap.read_dataset(ptype, "Velocity")[mask,:]*np.sqrt(self.a_0)  - bulkvel
            BH_dict['Mass'] = snap.read_dataset(ptype, "BH_Mass")[mask] /self.h
            BH_dict['Density'] = snap.read_dataset(ptype, "BH_Density")[mask] * self.h**2/self.a_0**3
            BH_dict['ParticleIDs'] = snap.read_dataset(ptype, "ParticleIDs")[mask]
            BH_dict['GroupNumber'] = snap.read_dataset(ptype, "GroupNumber")[mask]
            BH_dict['SmoothingLength'] = snap.read_dataset(ptype, "SmoothingLength")[mask] * self.a_0/self.h

        volume = {}
        volume['GroupNumber'] = groupnumber
        volume['CentreOfPotential'] = centre
        volume['BulkVelocity'] = bulkvel
        volume['r200'] = r200

        # Put everything together into one dictionary to return

        halo_dict = {}
        halo_dict['Gas'] = gas_dict
        halo_dict['DarkMatter'] = DM_dict
        halo_dict['Stars'] = star_dict
        if self.incBH:
            halo_dict['BlackHoles'] = BH_dict
        halo_dict['Volume'] = volume

        return halo_dict




    def save(self, halo_dict, halotype, path='/data6/arijdav1/saved_regions/', overwrite=False):  ###### For saving each subhalo as a seperate file

        # Change this line if you want to change the directory structure under which the results will be stored
        vol_dir = path + self.sim + '_' + self.run + '/' + self.tag + '/' + halotype + '/'

        # Uncomment this line if you want a more EAGLE-style directory structure without the 'halotype' descriptor
        #vol_dir = path + self.sim + '/' + self.run + '/' + self.tag + '/'

        vol_file = vol_dir + 'group' + str(halo_dict['Volume']['GroupNumber'])
        exist = os.path.exists(vol_file)

        # This is to prevent overwriting of existing data if you didn't check whether the files already exist
        if exist == True and overwrite == False:
            print 'These halos have already been saved. To overwrite, enter "overwrite" as a command line argument"'
            exit()

        # Generate the output file
        if exist == False or overwrite == True:
            if not os.path.exists(vol_dir):
                os.makedirs(vol_dir)
            f = h5py.File(vol_file + '.hdf5', 'w')
            vol = f.create_group('Volume')
            gas = f.create_group('Gas')
            dm = f.create_group('DarkMatter')
            stars = f.create_group('Stars')

            vol_dict = halo_dict['Volume']
            for k, key in enumerate(vol_dict.keys()):
                vol.create_dataset(key,data=vol_dict[key])

            gas_dict = halo_dict['Gas']
            for k, key in enumerate(gas_dict.keys()):
                gas.create_dataset(key,data=gas_dict[key])

            star_dict = halo_dict['Stars']
            for k, key in enumerate(star_dict.keys()):
                stars.create_dataset(key,data=star_dict[key])

            DM_dict = halo_dict['DarkMatter']
            for k, key in enumerate(DM_dict.keys()):
                dm.create_dataset(key,data=DM_dict[key])

            if self.incBH:
                bh = f.create_group('BlackHoles')
                BH_dict = halo_dict['BlackHoles']
                for k, key in enumerate(BH_dict.keys()):
                    bh.create_dataset(key,data=BH_dict[key])

            f.close()
