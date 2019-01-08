import eagle as E
import numpy as np
#import read as eread
import h5py
import os
from sys import exit
from tqdm import tqdm
from time import time
import numexpr as ne

class Snapshot_Volume(object):
    def __init__(self, centers, side_sizes, bulkvels, central_gnums,sgnums,
                 align_radius=0.03, 
                 run='L0100N1504', 
                 model='REFERENCE', 
                 tag='028_z000p000',
                 part_types = ['0','4'],
                 pdata_type = 'SNAPSHOT',
                 overwrite = True,
                 path = '/scratch/eagle/sav/',
                 halotype = 'any_halos',
                 timing_flag = False):
        vol_file = path+run+'_'+model+'/'+tag+'/FOF'+str(int(central_gnums[0]))+'_volume.hdf5'

        if overwrite == True or not os.path.exists(vol_file):    
            sim = '/data5/simulations/EAGLE/'+run+'/'+model+'/data'
            #get boxsize and particle positions
            boxsize = E.readAttribute(pdata_type, sim, tag, "/Header/BoxSize")
            h = E.readAttribute(pdata_type, sim, tag, "/Header/HubbleParam")
            Omega0 = E.readAttribute(pdata_type, sim, tag, "/Header/Omega0")
            OmegaLambda = E.readAttribute(pdata_type, sim, tag, "/Header/OmegaLambda")
            OmegaBaryon = E.readAttribute(pdata_type, sim, tag, "/Header/OmegaBaryon")
            a_0 = E.readAttribute(pdata_type, sim, tag, "/Header/ExpansionFactor")
            boxsize = boxsize/h
            masstable = E.readAttribute(pdata_type, sim, tag, "/Header/MassTable") / h

            self.mask_type = []
            self.pos_type = []
            self.vel_type = []
            self.mass_type = []
            self.density_type = []
            self.PID_type = []
            self.abund_type = []
            self.smooths = []
            self.gnum_type = []
            self.sfr_type = []
            self.metal = []
            self.temp = []

            num_subs = E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/NumOfSubhalos")
            r200 = E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/Group_R_Crit200")[num_subs>0]
            gns = central_gnums - 1
            self.r200 = r200[gns]

            #re-center to 'center' and remove particles outside box
            for t in part_types:


                masks = []
                pos = load_array('Coordinates', t, array_type=pdata_type, run=run, model=model, tag=tag)
                grouppos = []

                # Create a list of COPs where each COP is given in the co-ordinate system of the previous COP.
                shift_centres = np.zeros(np.shape(centers))
                shift_centres[0,:] = centers[0,:]
                for c in range(len(centers)-1):
                    shift_centres[c+1] = centers[c+1]-centers[c]

                inds = np.arange(len(pos[:,0]))
                npart = len(inds)

                print 'Creating apertures...'
                for ii in tqdm(range(len(shift_centres))):
                    centre = shift_centres[ii]

                    if ii == 0:
                        current_centre = np.array([0.,0.,0.])
                    else:
                        current_centre = centers[ii-1] # the current centre in the true box co-ordinate system
                    
                    #if (np.absolute(centre)+side_sizes[ii]/2.).any() > boxsize/2.: 
                    
                    if ((boxsize/2.+np.absolute(current_centre)) - (np.absolute(centre)+side_sizes[ii]/2.)).any() < 0.: # Is the group actually on the edge?
                        pos = ne.evaluate("pos-centre")
                        pos[pos[:,0]<(-1.*boxsize/2.),0] += boxsize
                        pos[pos[:,1]<(-1.*boxsize/2.),1] += boxsize
                        pos[pos[:,2]<(-1.*boxsize/2.),2] += boxsize


                        '''
                        pos -= (centre - boxsize/2.)
                        pos = np.mod(pos,boxsize)
                        pos -= boxsize/2.
                        '''
                        #if timing_flag:
                        print 'Wrapped box'

                    else: # Don't bother doing the wrapping if it doesn't affect the current group
                        s = time()
                        pos = ne.evaluate("pos-centre")
                        transform_time = time() - s
                    
                    rmax = side_sizes[ii]/2.

                    s = time()
                    # Chop along x
                    xs = pos[:,0]
                    mask = ne.evaluate('where(abs(xs)<rmax,True,False)').nonzero()[0]
                    cut_pos = pos[mask,:]
                    cut_inds = inds[mask]
                    # Along y
                    ys=cut_pos[:,1]
                    mask = ne.evaluate('where(abs(ys)<rmax,True,False)').nonzero()[0]
                    cut_pos = cut_pos[mask,:]
                    cut_inds = cut_inds[mask]
                    # Along z
                    zs = cut_pos[:,2]
                    mask = ne.evaluate('where(abs(zs)<rmax,True,False)').nonzero()[0]
                    cut_pos = cut_pos[mask,:]
                    cut_inds = cut_inds[mask]
                    chop_time = time() - s

                    s = time()
                    r2 = np.einsum('...j,...j->...',cut_pos,cut_pos) # get the radii from the centre
                    radial_calculation = time() - s

                    s = time()
                    mask = np.where(r2<rmax**2)[0] # make the mask
                    radial_masking = time() - s

                    masks.append(cut_inds[mask])
                    grouppos.append(cut_pos[mask])

                    if timing_flag:
                        print 'Time to transform box: ',transform_time
                        print 'Time to chop box up: ',chop_time
                        print 'Time to calculate radial distances: ',radial_calculation
                        print 'Time to mask by radius: ',radial_masking


                self.pos_type.append(grouppos)
                del pos
                del inds
                
                velarr = load_array('Velocity', t, array_type=pdata_type, run=run, model=model, tag=tag)
                groupvel = []
                for ii,mask in enumerate(masks):
                    vel = velarr[mask]
                    vel -= bulkvels[ii]
                    groupvel.append(vel)
                self.vel_type.append(groupvel)
                del velarr

                groupmass = []
                if t != '1':
                    massarr = load_array('Mass', t, array_type=pdata_type, run=run, model=model, tag=tag)
                if t == '1':
                    massarr = np.ones(npart)*masstable[1]
                for ii,mask in enumerate(masks):
                    groupmass.append(massarr[mask])
                self.mass_type.append(groupmass)
                del massarr

                if t in ['0','4','5']:
                    groupdensity = []
                    if t == '0':
                        densityarr = load_array('Density', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    if t == '4':
                        densityarr = load_array('BirthDensity', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    if t == '5':
                        densityarr = load_array('BH_Density', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    for ii,mask in enumerate(masks):
                        groupdensity.append(densityarr[mask])
                    self.density_type.append(groupdensity)
                    del densityarr

                grouppids = []
                PIDs = load_array('ParticleIDs', t, array_type=pdata_type, run=run, model=model, tag=tag)
                for ii,mask in enumerate(masks):
                    grouppids.append(PIDs[mask])
                self.PID_type.append(grouppids)
                del PIDs

                groupgnums = []
                gnums = load_array('GroupNumber', t, array_type=pdata_type, run=run, model=model, tag=tag)
                for ii,mask in enumerate(masks):
                    groupgnums.append(gnums[mask])
                self.gnum_type.append(groupgnums)
                del gnums
                
                if t == '0':
                    groupsfrs = []
                    sfr = load_array('StarFormationRate', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    for ii,mask in enumerate(masks):
                        groupsfrs.append(sfr[mask])
                    self.sfr_type = groupsfrs

                    grouptemps = []
                    temps = load_array('Temperature', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    for ii,mask in enumerate(masks):
                        grouptemps.append(temps[mask])
                    self.temp = grouptemps

                if t == '4':
                    groupsftimes = []
                    sftimes = load_array('StellarFormationTime', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    for ii,mask in enumerate(masks):
                        groupsftimes.append(sftimes[mask])
                    self.sftimes = groupsftimes

                if t in ['0','4']:
                    H = load_array('SmoothedElementAbundance/Hydrogen', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    He = load_array('SmoothedElementAbundance/Helium', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    C = load_array('SmoothedElementAbundance/Carbon', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    N = load_array('SmoothedElementAbundance/Nitrogen', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    O = load_array('SmoothedElementAbundance/Oxygen', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    Ne = load_array('SmoothedElementAbundance/Neon', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    Mg = load_array('SmoothedElementAbundance/Magnesium', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    Si = load_array('SmoothedElementAbundance/Silicon', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    #S = load_array('SmoothedElementAbundance/Sulphur', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    #Ca = load_array('SmoothedElementAbundance/Calcium', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    Fe = load_array('SmoothedElementAbundance/Iron', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    abunds = np.dstack((H,He,C,N,O,Ne,Mg,Si,Fe))[0]
                    del H,He,C,N,O,Ne,Mg,Si,Fe
                    groupabunds = []
                    for ii,mask in enumerate(masks):
                        groupabunds.append(abunds[mask])
                    self.abund_type.append(groupabunds)
                    del abunds

                    groupmetal = []
                    metal = load_array('SmoothedMetallicity', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    for ii,mask in enumerate(masks):
                        groupmetal.append(metal[mask])
                    self.metal.append(groupmetal)
                    del metal


                if t in ['0','4','5']:
                    groupsmooths = []
                    smooths = load_array('SmoothingLength', t, array_type=pdata_type, run=run, model=model, tag=tag)
                    for ii,mask in enumerate(masks):
                        groupsmooths.append(smooths[mask])
                    self.smooths.append(groupsmooths)
                    del smooths
            self.mask_type.append(masks)
                
    	    for n in range(len(part_types)):
                if part_types[n] == '4':
                    print 'Stars are part index ',n
                    si = n
                else:
                    continue
                self.transform = []
                for ii in range(0,len(self.pos_type[0])):
                    #perform alignment to Jz of star particles in a given R
                    Rs = np.sqrt(self.pos_type[si][ii][:,0]**2+self.pos_type[si][ii][:,1]**2+self.pos_type[si][ii][:,2]**2)
                    radmask = Rs < align_radius
                    starj = np.cross(self.pos_type[si][ii],self.vel_type[si][ii]*self.mass_type[si][ii][:,np.newaxis])
                    r200j = starj[radmask]
                    tot_ang_mom = np.sum(r200j, axis = 0)
                    a = np.matrix([tot_ang_mom[0],tot_ang_mom[1],tot_ang_mom[2]])/np.linalg.norm([tot_ang_mom[0],tot_ang_mom[1],tot_ang_mom[2]])
                    b = np.matrix([0,0,1])
                    v = np.cross(a,b)
                    s = np.linalg.norm(v)
                    c = np.dot(a,b.T)
                    vx = np.matrix([[0,-v[0,2],v[0,1]],[v[0,2],0,-v[0,0]],[-v[0,1],v[0,0],0]])
                    transform = np.eye(3,3) + vx + (vx*vx)*(1/(1+c[0,0]))
                    self.transform.append(transform)
            self.run = run
            self.model = model
            self.tag = tag
            self.halotype = halotype
            self.centers = centers
            self.side_sizes = side_sizes
            self.bulkvels = bulkvels
            self.central_gnums = central_gnums
    	    self.part_types = part_types

	    
            print 'Applying transformations...'
            for ii,ptype in enumerate(part_types):
                for jj,transform in tqdm(enumerate(self.transform)):
                    try:
                        self.pos_type[ii][jj] = np.array([np.dot(transform,self.pos_type[ii][jj][i].T) for i in range(0,len(self.pos_type[ii][jj]))])[:,0]
                        self.vel_type[ii][jj] = np.array([np.dot(transform,self.vel_type[ii][jj][i].T) for i in range(0,len(self.vel_type[ii][jj]))])[:,0]
                    except IndexError:
                        continue

    def save_individuals(self, path='/data5/arijdav1/saved_regions/', overwrite=False, addDM=False):      ###### For saving each subhalo as a seperate file
        for jj,center in enumerate(self.centers):
            vol_file = path+self.run+'_'+self.model+'/'+self.tag+'/'+self.halotype+'/group'+str(self.central_gnums[jj])+'.hdf5'
            exist = os.path.exists(vol_file)

            if exist == True and overwrite == False and not addDM:
                print 'These halos have already been saved. To overwrite, enter "overwrite" as a command line argument"'
                exit()

            if exist == True and overwrite == False and addDM:
                f = h5py.File(vol_file, 'r+')
                # as defined in run_snaptools, the part_type index for the DM is 0 as we only retrive DM and stars for 'addDM'
                try:
                    f['DarkMatter'].create_dataset('Coordinates',data=self.pos_type[0][jj])
                    f['DarkMatter'].create_dataset('Velocity',data=self.vel_type[0][jj])
                    f['DarkMatter'].create_dataset('Mass',data=self.mass_type[0][jj])
                    f['DarkMatter'].create_dataset('PID',data=self.PID_type[0][jj])
                    f['DarkMatter'].create_dataset('GroupNumber',data=self.gnum_type[0][jj])
                    f.close()
                except RuntimeError:
                    f.close()
                    continue
                

            if exist == False or overwrite==True:
                ensure_dir(vol_file)
                print 'Saving...'
                f = h5py.File(vol_file, 'w')
                vol = f.create_group('Volume')
                gas = f.create_group('Gas')
                dm = f.create_group('DarkMatter')
                stars = f.create_group('Stars')
                bh = f.create_group('BlackHoles')
                types = []
                if '0' in self.part_types:
                    types.append(gas)
                if '1' in self.part_types:
                    types.append(dm)
                if '4' in self.part_types:
                    types.append(stars)
                if '5' in self.part_types:
                    types.append(bh)
                aind = 0
                sind = 0
                for ii,t in enumerate(types):
                    t.create_dataset('Coordinates', data=self.pos_type[ii][jj])
                    t.create_dataset('Velocity', data=self.vel_type[ii][jj])
                    t.create_dataset('Mass', data=self.mass_type[ii][jj])
                    t.create_dataset('PID', data=self.PID_type[ii][jj])
                    t.create_dataset('GroupNumber', data=self.gnum_type[ii][jj])
                    
                    if self.part_types[ii] in ['0']:
                        t.create_dataset('StarFormationRate', data=self.sfr_type[jj])
                        t.create_dataset('Temperature', data=self.temp[jj])

                    if self.part_types[ii] in ['4']:
                        t.create_dataset('StellarFormationTime', data=self.sftimes[jj])
                    
                    if self.part_types[ii] in ['0','4']:
                        t.create_dataset('Abundances', data=self.abund_type[aind][jj])
                        t.create_dataset('Metallicity', data=self.metal[aind][jj])
                        aind += 1

                    if self.part_types[ii] in ['0','4','5']:
                        t.create_dataset('SmoothingLength', data=self.smooths[sind][jj])
                        t.create_dataset('Density', data=self.density_type[sind][jj])
                        sind += 1
                
                vol.create_dataset('r200', data=self.r200[jj])
                vol.create_dataset('Run', data=self.run)
                vol.create_dataset('Model', data=self.model)
                vol.create_dataset('Tag', data=self.tag)
                vol.create_dataset('Centers', data=self.centers[jj])
                vol.create_dataset('SideSize', data=self.side_sizes[jj])
                vol.create_dataset('BulkVelocities', data=self.bulkvels[jj])
                vol.create_dataset('CentralGroupNumbers', data=self.central_gnums[jj])
                vol.create_dataset('Transforms', data=self.transform[jj])
                f.close()


class fofinfo(object):
    def __init__(self, gnum,  run='L0100N1504', model='REFERENCE', tag='028_z000p000'):
        self.run, self.model, self.tag = run, model, tag
        sim = '/data5/simulations/EAGLE/'+run+'/'+model+'/data'
        fsid = np.array(E.readArray("SUBFIND_GROUP" , sim, tag, "/FOF/FirstSubhaloID"))
        groupnumber = np.array(E.readArray("SUBFIND" , sim, tag, "/Subhalo/GroupNumber"))[fsid]
        subgroupnumber = np.array(E.readArray("SUBFIND", sim, tag, "/Subhalo/SubGroupNumber"))[fsid]
        self.gnum = gnum
        ind = np.where(np.in1d(groupnumber, gnum) & (subgroupnumber == 0))
        print ind
        self.CoP = np.array(E.readArray("SUBFIND", sim, tag, "/Subhalo/CentreOfPotential"))[fsid][ind]
        self.subhalovel = np.array(E.readArray("SUBFIND", sim, tag, "/Subhalo/Velocity"))[fsid][ind]
        self.r_200 = np.array(E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/Group_R_Crit200"))[ind]
        self.m_200 = np.array(E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/Group_M_Crit200")*1e10)[ind]
        self.tot_ang_mom = np.array(E.readArray("SUBFIND", sim, tag, "/Subhalo/Stars/Spin"))[fsid][ind]
        self.stellar_mass = np.array(E.readArray("SUBFIND", sim, tag, "/Subhalo/Stars/Mass") * 1e10)[fsid][ind]
    def save_fofinfo(self, path = '/gal/eagle/sav/', overwrite = True):
        fof_file = path+self.run+'_'+self.model+'/'+self.tag+'/FOF'+str(int(self.gnum))+'_fof.hdf5'
        ensure_dir(fof_file)
        exist = os.path.exists(fof_file)
        if exist == False or overwrite==True:
            print 'saving...'
            f = h5py.File(fof_file, 'w')
            vol = f.create_group('FoF')
            vol.create_dataset('Run', data=self.run)
            vol.create_dataset('Model', data=self.model)
            vol.create_dataset('Tag', data=self.tag)
            vol.create_dataset('GroupNumber', data=self.gnum)
            vol.create_dataset('CoP', data=self.CoP)
            vol.create_dataset('BulkVelocity', data=self.subhalovel)
            vol.create_dataset('R_200', data=self.r_200)
            vol.create_dataset('M_200', data=self.m_200)
            vol.create_dataset('Spin', data=self.tot_ang_mom)
            vol.create_dataset('StellarMass', data=self.stellar_mass)
            f.close()
'''
def load_info(run,model,tag, arr_type='SNAP'):
    return eread.hdf5_info(run,model,tag)
'''
def load_array(array_label, parttype, array_type='SNAP', run='L0100N1504', model='REFERENCE', tag='028_z000p000'):
    sim = '/data5/simulations/EAGLE/'+run+'/'+model+'/data'
    return np.array(E.readArray(array_type, sim, tag, '/PartType'+parttype+'/'+array_label))
    
def load_attribute(att_label, array_type='SNAP', run='L0100N1504', model='REFERENCE', tag='028_z000p000'):
    sim = '/data5/simulations/EAGLE/'+run+'/'+model+'/data'
    return E.readAttribute(att_label, sim, tag, '/Header/'+att_label)
    

def abundratios(abunds):
    H, O, Mg, Si, Fe = abunds[:,0], abunds[:,1], abunds[:,2], abunds[:,3], abunds[:,4]
    solar_H = 0.706498
    solar_Fe = 0.00110322
    solar_O = 0.00549262
    solar_Si = 0.000682587
    solar_Mg = 0.000590706
    O_H, Mg_H, Si_H, Fe_H = np.log10(O/H)-np.log10(solar_O/solar_H), np.log10(Mg/H)-np.log10(solar_Mg/solar_H), np.log10(Si/H)-np.log10(solar_Si/solar_H), np.log10(Fe/H)-np.log10(solar_Fe/solar_H)
    a_Fe = (O_H+Mg_H+Si_H)/3-Fe_H
    return Fe_H, a_Fe
    
def ensure_dir(f):
    """ Ensure a a file exists and if not make the relevant path """
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
