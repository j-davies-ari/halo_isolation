import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def pload(fname): # Loads the pickle file for a particular data sample, given a filename
    infile = '/home/arijdav1/Desktop/query_results/%s.pkl'%(fname)
    pkl_file = open(infile, 'rb')
    myData = pickle.load(pkl_file)
    pkl_file.close()
    
    return myData
    
def tmload(fname): # Loads the pickle file for a particular data sample, given a filename
    infile = '/home/arijdav1/Desktop/query_results/%s.pkl'%(fname)
    pkl_file = open(infile, 'rb')
    myData = pickle.load(pkl_file)
    pkl_file.close()
    
    return myData
    
def qprint(fname): # Prints the query used to generate a data sample, given the filename of that sample
    infile = '/home/arijdav1/Desktop/query_results/%s_query.pkl'%(fname)
    query = open(infile,'r')
    #print query.read()
    for line in query:
        print line,
    
def prop(dic): # Prints the galaxy properties available in the data sample
    for i in dic:
        print i,', '
        
def zlist(): # Gives a list of the snapshot redshifts in EAGLE up to z~10
    return np.loadtxt('/Users/jonathan/Desktop/Masters_Project/Data/z_list.dat')  
    
def get_binedges(z_store): # Creates an array of bin edges for variable-width histogram plotting, given the bin centres
    bins = []
    bins.append(z_store[0]-(z_store[1]-z_store[0])/2)
    for j in range(len(z_store)):
        if j+1 == len(z_store):
            break
        else:
            bins.append(z_store[j]+(z_store[j+1]-z_store[j])/2)
    bins.append(z_store[-1]+(z_store[-1]-z_store[-2])/2)
    
    return bins
    
def get_bincentres(binedges): # Finds the centre points of a set of bin edges
    bincentres = []
    for i in range(len(binedges)):
        if i+1 == len(binedges):
            break
        else:
            bincentres.append((binedges[i+1]+binedges[i])/2)
    return bincentres
    
def get_binsizes(binedges):
    binsizes = []
    for i in range(len(binedges)):
        if i+1 == len(binedges):
            break
        else:
            binsizes.append(binedges[i+1]-binedges[i])
    return binsizes
    
def dmass(mass,metal): # Returns the dust masses for a galaxy sample
    return 0.4*metal*mass
    
def submmflux(sfr,mass,metal):
    ngal = len(sfr)
    ones = np.ones(ngal)
    fl = 0.81*(sfr/(ones*100))**0.43*(dmass(mass,metal)/(ones*10**8))**0.58
    return fl
    
def submm_2013_flux(sfr,mass,metal):
    ngal = len(sfr)
    ones = np.ones(ngal)
    fl = 0.81*(sfr/(ones*100))**0.43*(dmass(mass,metal)/(ones*10**8))**0.54
    return fl
    
def pop_meanplot(sfrm,subm,sfrrefm,subrefm):
    plt.axvline(x=sfrm,linewidth=1,color='red',ls=':')
    plt.axvline(x=subm,linewidth=1,color='blue',ls=':')
    plt.axvline(x=sfrrefm,linewidth=1,color='magenta',ls=':')
    plt.axvline(x=subrefm,linewidth=1,color='cyan',ls=':')
    
def find_sigma(array,upper,lower): # Finds indices at which one standard
    up = (np.abs(array-upper)).argmin() # deviation occurs, given cumulative data
    dn = (np.abs(array-lower)).argmin()
    return up, dn
    
def find_perc(a,percs): # Finds intervals at which certain percentiles are
    inds = []            # first reached, fiven cumulative data and an array
    for p in percs:       # of percentiles of interest
        tup = np.where(a>=p)
        locs = tup[0]
        if locs.size == 0:
            break
        inds.append(locs[0]) 
    return inds

def find_xval(array,xval): # Finds closest index to which a given x value occurs
        return np.abs(array-xval).argmin()

def z_thirds(z): # Splits redshift data into 3 by number and returns the
    for p in range(1000):
        if z[-1]-z[-2] > 1:
            tmp = np.around(z[-1],decimals=2)
            z = np.delete(z,-1)
            print 'Outlier galaxy removed at z = ',tmp
        elif z[-1]<= z[-2]:
            break
    tot = len(z) # redshift ranges of the first and third thirds.
    thd = np.floor(np.float(tot)/3.)
    first = z[0:thd]
    second = z[thd:(2*thd)]
    third = z[(2*thd):]  
    if len(third) > len(second):
        second = np.append(second,[third[0]])
        third = np.delete(third,[0])
    lowz = np.array([first[0],first[-1]])
    highz = np.array([third[0],third[-1]])
    return lowz, highz
    
def prob_stdev(n_sel,low,high):
    p_low = low**(1/n_sel)
    p_high = high**(1/n_sel)
    return p_low, p_high
    
def prob_lines(ratio,n_sel):
    return ratio**(n_sel)

def plot_scatter(x,y,shade,xlim,ylim,xlabel,ylabel,clabel,xls,yls,cls,save,set_bounds):
        plt.figure(figsize=(10,8))
        col = cm.get_cmap('Reds')
        colr = cm.get_cmap('Blues')
        plt.scatter(x[1],y[1],c=shade[1],cmap=colr,marker='o',s=15)
        plt.scatter(x[0],y[0],s=60,c=shade[0],marker='s',cmap=col)
        plt.xlabel(xlabel,fontsize=xls)
        plt.ylabel(ylabel,fontsize=yls)
        if set_bounds:
            plt.xlim(xlim[0],xlim[1])
            plt.ylim(ylim[0],ylim[1])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(clabel,fontsize=cls)
        plt.savefig(save)
        plt.show()
        
def scatter_highlight_desc(x,y,shade,desc_inds,hsf_ids,xlim,ylim,xlabel,ylabel,clabel,xls,yls,cls,save,set_bounds,show_ids):
        plt.figure(figsize=(10,8))
        col = cm.get_cmap('Reds')
        colr = cm.get_cmap('Blues')
        xs = x[0]
        ys = y[0]
        rxs = x[1]
        rys = y[1]
        plt.scatter(rxs,rys,c=shade[1],cmap=colr,marker='o',s=15)
        plt.scatter(xs,ys,s=40,c=shade[0],cmap=col)

        if show_ids:        
            for j in range(len(hsf_ids)):
                plt.text(xs[j],ys[j],'%s'%(str(hsf_ids[j])))
        
        plt.xlabel(xlabel,fontsize=xls)
        plt.ylabel(ylabel,fontsize=yls)
        if set_bounds:
            plt.xlim(xlim[0],xlim[1])
            plt.ylim(ylim[0],ylim[1])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(clabel,fontsize=cls)
        plt.scatter(xs[desc_inds],ys[desc_inds],edgecolor='g',facecolor='none',s=160,marker='o')
        plt.savefig(save)
        plt.show()        
        
def scatter_follow(x,y,shade,desc_inds,hsf_ids,foll,xlim,ylim,xlabel,ylabel,clabel,xls,yls,cls,save,set_bounds,show_ids):
        plt.figure(figsize=(10,8))
        col = cm.get_cmap('Reds')
        colr = cm.get_cmap('Blues')
        xs = x[0]
        ys = y[0]
        rxs = x[1]
        rys = y[1]
        plt.scatter(rxs,rys,c=shade[1],cmap=colr,marker='o',s=15)
        plt.scatter(xs,ys,s=40,c=shade[0],cmap=col)
        
        if show_ids:
            for j in range(len(hsf_ids)):
                plt.text(xs[j],ys[j],'%s'%(str(hsf_ids[j])))
        
        plt.xlabel(xlabel,fontsize=xls)
        plt.ylabel(ylabel,fontsize=yls)
        if set_bounds:
            plt.xlim(xlim[0],xlim[1])
            plt.ylim(ylim[0],ylim[1])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(clabel,fontsize=cls)
        plt.scatter(xs[desc_inds],ys[desc_inds],edgecolor='g',facecolor='none',s=160,marker='o')
        plt.scatter(xs[foll],ys[foll],edgecolor='magenta',facecolor='none',s=200,marker='D')
        plt.savefig(save)
        plt.show()    

        
def snap_select(z,index):
    snap_zs = zlist()
    zsel = snap_zs[-index]
    lz = zsel-0.001
    hz = zsel+0.001
    iz = np.where((z>=lz) & (z<=hz))
    iz_out = iz[0][::-1]
    return np.around(zsel,decimals=2), iz_out

def ourstat(sample_xvals,ratio,cont_xvals):
        N = len(sample_xvals)
        single_results = []
        for i in range(N):
            frac = ratio[find_xval(cont_xvals,sample_xvals[i])]
            single_results.append(frac*(1-frac))
        #plt.figure()
        #plt.scatter(sample_xvals,single_results)
        #plt.xlabel(r'Our statistic')
        #plt.show()
        return (1/np.float(N))*np.sum(single_results)
        
def shannon(sample_xvals,ratio,cont_xvals):
        N = len(sample_xvals)
        single_results = []
        for i in range(N):
            frac = ratio[find_xval(cont_xvals,sample_xvals[i])]
            single_results.append(frac*np.log10(frac))
        #plt.figure()
        #plt.scatter(sample_xvals,single_results)
        #plt.show()
        #plt.xlabel(r'Shannon Entropy')
        return -(1/np.float(N))*np.sum(single_results)
        
        
def link_tree(xaxis,yaxis,z,gal,des): # Links a meger tree, specific to mergertree.py
    for i in range(len(z)):
        if gal[i] == 21268427:
            continue
        desc_idx = np.where(gal==des[i])
        dx = desc_idx[0]
        #print dx
        #if 1 < gm[i]/gm[dx] < 4:
        #    print 'Major merger!'
        #    plt.scatter(xaxis[i],yaxis[i],edgecolor='m',facecolor='none',s=160,marker='o')
        #    plt.scatter(xaxis[dx],yaxis[dx],edgecolor='m',facecolor='none',s=160,marker='o')
        if dx.size == 0:
            desc_idx = np.where(gal==gal[i])
            dx = desc_idx[0]
        #    print dx,' nomerge'
            #continue
        if dx.size==0:
            continue
        xarr = np.array([xaxis[i],xaxis[dx]])
        yarr = np.array([yaxis[i],yaxis[dx]])
        #print xarr
        #print yarr
        plt.plot(xarr,yarr,color='black')

def link_mainbranch(gal_z,prog_z,gal_id,prog_id): # Links a main-branch and merger tree
    next_z = 0
    for i in range(len(gal_z)):
        old_z = next_z
        for j in range(100000):
            if prog_z[j]>gal_z[i]:
                next_z =prog_z[j]
                break
        if next_z-old_z > 0.09:
            continue
        print next_z
        li = np.where(prog_z==next_z)
        links = li[0]
        startzs = np.ones(len(links))*gal_z[i]
        startids = np.ones(len(links))*gal_id[i]
        endzs = prog_z[links]
        endids = prog_id[links]
        for k in range(len(links)):
            xarr = np.array([startids[k],endids[k]])
            yarr = np.array([startzs[k],endzs[k]])
            plt.plot(xarr,yarr,c='black',zorder=1)
    plt.plot(np.array([gal_id[0],gal_id[-1]]),np.array([gal_z[0],gal_z[-1]]),zorder=1,c='k')
    
    
def ax_link_mainbranch(ax,gal_z,prog_z,gal_id,prog_id): # Links a main-branch and merger tree
    next_z = 0
    for i in range(len(gal_z)):
        old_z = next_z
        for j in range(100000):
            if prog_z[j]>gal_z[i]:
                next_z =prog_z[j]
                break
        if next_z-old_z > 0.09:
            continue
        print next_z
        li = np.where(prog_z==next_z)
        links = li[0]
        startzs = np.ones(len(links))*gal_z[i]
        startids = np.ones(len(links))*gal_id[i]
        endzs = prog_z[links]
        endids = prog_id[links]
        for k in range(len(links)):
            xarr = np.array([startids[k],endids[k]])
            yarr = np.array([startzs[k],endzs[k]])
            ax.plot(xarr,yarr,c='black',zorder=1)
    ax.plot(np.array([gal_id[0],gal_id[-1]]),np.array([gal_z[0],gal_z[-1]]),zorder=1,c='k')





'''
# BACKUP    
def link_mainbranch(gal_z,prog_z,gal_id,prog_id): # Links a main-branch and merger tree
    for i in range(len(gal_z)):
        xarr = np.array([gal_id[i],prog_id[i]])
        yarr = np.array([gal_z[i],prog_z[i]])
        plt.plot(xarr,yarr,c='black')        
'''
    
def unique_entries(arr): # Returns the indices of the first instance of a unique element in a 1D array
    output = []
    val = arr[0]
    output.append(0)
    for i in range(len(arr)):
        if arr[i] != val:
            val = arr[i]
            output.append(i)
    return output
    
def standard_plane(colouring,clabel,outfile):
    figure, ax1 = plt.subplots(1,figsize=(10,8))
    col = cm.get_cmap('Reds')
    colr = cm.get_cmap('Greens')
    refsc = ax1.scatter(all_xs,all_ys,c=colouring,cmap=colr,marker='s',edgecolor = 'none',s=160)
    scat = ax1.scatter(hsf_lagos_sfr_zpredict,hsf_lagos_sfr_eagle,c=z,cmap=col,marker='o',s=60)
    plt.plot(np.arange(-5,5),np.arange(-5,5),c='grey')    
    ignore = ax1.scatter(np.arange(-100,-100.0001),np.arange(-100,-100.0001),c='r',marker='s',s=60,label=r'$\mathrm{SFR}$ $>$ $60$ $M_{\odot}$ $\mathrm{yr}^{-1}$, $%s<z<%s$'%(str(lzr),str(hzr))) 
    ignore2 = ax1.scatter(np.arange(-100,-100.0001),np.arange(-100,-100.0001),c='g',marker='s',s=60,label=r'$\mathrm{All}$, $%s<z<%s$'%(str(lzr),str(hzr))) 
    plt.plot(bins,predict_med,c='b',lw=2)
    plt.plot(bins,predict_med+std,c='b',lw=2,ls='--')
    plt.plot(bins,predict_med-std,c='b',lw=2,ls='--')    
    ax1.set_xlim(-4,3)
    ax1.set_ylim(-4,3)
    if k == 0 and sm_cut:
        ax1.set_xlim(-1.5,2.5)
        ax1.set_ylim(-1.5,2.5)
    ax1.set_ylabel(ylab, fontsize=18)
    ax1.set_xlabel(xlab, fontsize=18)
    ax1.legend(frameon=False,loc=2,prop={'size':12})
    cbar1 = figure.colorbar(refsc,ax=ax1,shrink=0.95)
    cbar1.ax.set_ylabel(clabel,fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0,hspace=0)
    if k == 0:
        which = 'highz'
    elif k == 1:
        which = 'lowz'
    plt.savefig(fname)
    plt.show()
    
