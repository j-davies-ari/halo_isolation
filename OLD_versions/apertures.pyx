def create_apertures(posarr,centers,side_sizes,boxsize):
    masks = []
    grouppos = []
    print 'Creating apertures...'
    for ii,center in enumerate(centers):
        maxrad = side_sizes[ii]/2.

        if (abs(center)+maxrad).any() > boxsize/2.: # Is the group actually on the edge?
            pos = posarr - (center - boxsize/2.)
            pos %= boxsize
            pos -= boxsize/2.

        else: # Don't bother doing the wrapping if it doesn't affect the current group
            pos = posarr - center

        r = (pos[:,0]**2+pos[:,1]**2+pos[:,2]**2)**0.5

        mask = []
        for i in range(len(r)):
            if r[i] < maxrad:
                mask.append(i)

        #mask = np.where(r<maxrad)[0] # make the mask
        masks.append(mask)
        grouppos.append(pos[mask])
        print 'Done halo ',ii
    return grouppos, masks
