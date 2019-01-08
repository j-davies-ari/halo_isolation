import eagle as E
import numpy as np
import snaptools as snaptools
from sys import exit
from sys import argv
import os
import eagleSqlTools as sql
import pickle

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

snapdict = dict([('28',['028_z000p000',0.]),\
                ('27',['027_z000p101',0.101]),\
                ('26',['026_z000p183',0.183]),\
                ('25',['025_z000p271',0.271]),\
                ('24',['024_z000p366',0.366]),\
                ('23',['023_z000p503',0.503]),\
                ('22',['022_z000p615',0.615]),\
                ('21',['021_z000p736',0.736]),\
                ('20',['020_z000p865',0.865]),\
                ('19',['019_z001p004',1.004]),\
                ('18',['018_z001p259',1.259]),\
                ('17',['017_z001p487',1.487]),\
                ('16',['016_z001p737',1.737]),\
                ('15',['015_z002p012',2.012]),\
                ('14',['014_z002p237',2.237]),\
                ('13',['013_z002p478',2.478]),\
                ('12',['012_z003p017',3.017]),\
                ('11',['011_z003p528',3.528]),\
                ('10',['010_z003p984',3.984]),\
                ('9',['009_z004p485',4.485]),\
                ('8',['008_z005p037',5.037]),\
                ('7',['007_z005p487',5.487]),\
                ('6',['006_z005p971',5.971]),\
                ('5',['005_z007p050',7.05]),\
                ('4',['004_z008p075',8.075]),\
                ('3',['003_z008p988',8.988]),\
                ('2',['002_z009p993',9.993]),\
                ('1',['001_z015p132',15.132]),\
                ('0',['000_z020p000',20.])])

snapinfo = dict([('28',0.),\
                ('27',0.101),\
                ('26',0.183),\
                ('25',0.271),\
                ('24',0.366),\
                ('23',0.503),\
                ('22',0.615),\
                ('21',0.736),\
                ('20',0.865),\
                ('19',1.004),\
                ('18',1.259),\
                ('17',1.487),\
                ('16',1.737),\
                ('15',2.012),\
                ('14',2.237),\
                ('13',2.478),\
                ('12',3.017),\
                ('11',3.528),\
                ('10',3.984),\
                ('9',4.485),\
                ('8',5.037),\
                ('7',5.487),\
                ('6',5.971),\
                ('5',7.05),\
                ('4',8.075),\
                ('3',8.988),\
                ('2',9.993),\
                ('1',15.132),\
                ('0',20.)])

##################################################################
# Choose what you want
simulation = 'L0100N1504'
run = 'REFERENCE'
simlabel = 'RefL0100N1504'

#redshift_labels = ['0.0','0.5','1.0','2.0','5.0']
redshift_labels = ['2.0',]
#snapnums = [28,23,19,15,8]
snapnums = [15,]

#redshift_labels = ['0.0']
#redshift_labels = ['0.5','2.0','3.0']

galtype = 'all' # 'dwarf', 'Lstar', 'group' or 'all'

particle_types = ['0','4']

##################################################################

outfile = '/data5/arijdav1/merger_data/radial_profiling/mergertrees_'+simulation+'_'+run+'.pkl'

if not os.path.exists(outfile) or 'requery' in argv:
    # Look up which galaxies to get using the merger trees
    print 'Connecting to EAGLE database...'
    con = sql.connect("jdavies", password="L113Cg61")
    print 'Connected.'

    print 'Executing queries...'

    # Define mass bins, descending in mass so that we roughly ascend in group number
    mass_bins = np.log10(np.logspace(11.75,13.,60))[::-1]

    for m in range(len(mass_bins)-1):
        myQuery = '''-- Select the quantities we want
                    SELECT
                    DES.GroupNumber as GroupNumber,
                    PROG.SnapNum as prog_SnapNum,
                    PROG.GroupNumber as prog_GroupNumber,
                    FOF.Group_M_Crit200 as M200

                    -- Define aliases for the three tables
                    FROM
                    %s_Subhalo as PROG,
                    %s_Subhalo as DES,
                    %s_FOF as FOF

                    -- Apply the conditions
                    WHERE
                    FOF.Group_M_Crit200 between POWER(10.,%s) and POWER(10.,%s) -- Select halo mass
                    and FOF.GroupID = DES.GroupID
                    and DES.SubGroupNumber = 0
                    and DES.SnapNum = 28 -- At redshift z=0

                    and PROG.SnapNum IN ('28','23','19','15','8')
                    and PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID -- Then find galaxy progenitors

                    ORDER BY
                    DES.GroupNumber asc, PROG.Redshift asc'''%((simlabel,simlabel,simlabel,str(mass_bins[m+1]),str(mass_bins[m])))

        Data = sql.execute_query(con, myQuery)

        '''
        print Data['GroupNumber']
        print Data['M200']
        print Data['prog_SnapNum']
        print Data['prog_GroupNumber']
        '''
        
        try:
            gns = np.hstack((gns,np.array(Data['GroupNumber'])))
            M200s = np.hstack((M200s,np.array(Data['M200'])))
            prog_snapnums = np.hstack((prog_snapnums,np.array(Data['prog_SnapNum'])))
            prog_groupnums = np.hstack((prog_groupnums,np.array(Data['prog_GroupNumber'])))
            
        except NameError:
            gns = np.array(Data['GroupNumber'])
            M200s = np.array(Data['M200'])
            prog_snapnums = np.array(Data['prog_SnapNum'])
            prog_groupnums = np.array(Data['prog_GroupNumber'])

        #print gns

        print 'Bin ',m,' - ',len(np.unique(gns)),' galaxy main branches found'


    out_data = dict([('GroupNumber',np.asarray(gns)),\
                        ('final_M200',np.asarray(M200s)),\
                        ('prog_SnapNum',np.asarray(prog_snapnums)),\
                        ('prog_GroupNumber',np.asarray(prog_groupnums))])

    output = open(outfile, 'w')

    pickle.dump(out_data,output)

    output.close()

    exit()

else:
    infile = open(outfile,'rb')
    Data = pickle.load(infile)
    infile.close()
    gns = np.array(Data['GroupNumber'])
    M200s = np.array(Data['final_M200'])
    prog_snapnums = np.array(Data['prog_SnapNum'])
    prog_groupnums = np.array(Data['prog_GroupNumber'])


if addDM:
    particle_types = ['1','4']

for k in range(len(redshift_labels)):
    print 'Running halo isolator for ',simulation,'-',run
    print 'Redshift: ',redshift_labels[k]

    current_snap = snapnums[k]

    groupnums = prog_groupnums[prog_snapnums==current_snap]

    print groupnums


    snap = snapdict[str(current_snap)][0]
    
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

    r200s = r200s[groupnums]
    centres = subfind_centres[groupnums]
    subgroupnums = subfind_subgroupnums[groupnums]
    bulkvels= subfind_bulkvels[groupnums]

    print 'Scanning existing galaxy files...'
    nonexist_indices = []
    for c in range(len(centres)):
        savepath='/data5/arijdav1/mergertree_regions/'+simulation+'_'+run+'/'+snap+'/'+galtype+'/group'+str(groupnums[c])+'.hdf5'
        if not os.path.exists(savepath):
            nonexist_indices.append(c)

    if not overwrite and not addDM:
        if not nonexist_indices:
            print 'All selected galaxies have already been saved. To overwrite, enter "overwrite" as a command line argument'
            continue

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

    volume.save_individuals(path='/data5/arijdav1/mergertree_regions/',overwrite=overwrite,addDM=addDM)


