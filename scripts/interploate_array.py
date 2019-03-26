def interpolate(X,inflate):
    '''
    Generate a copy of the array, inflated in all dimensions by provided ratio
        currently works for a 3 dimensional array
    '''
    # determine shape of new array (-1 and +1 for the edges)
    shape = (np.array(normal.shape) - 1 ) * inflate +1
    # generate empty array
    inter = np.zeros(shape=shape)
    # insert values from the original array, equally spaced
    inter[tuple(slice(None,None,inflate) for i in range(X.ndim))] = X    
    # loop through each dimension, and calculate inflated indices
    for n_orig in range(X.shape[0]):
        n = n_orig * inflate
        for k_orig in range(X.shape[1]):
            k = k_orig*inflate
            for i_orig in range(X.shape[2]-1):
                i_new = i_orig*inflate
                j_new = (i_orig+1)*inflate
                # interpolate along the third dimension linearly
                inter[n,k][i_new:j_new+1] = np.linspace(inter[n,k][i_new],inter[n,k][j_new],inflate+1)
        # fill gaps between the interpolated values through the other dimensions
        for i in range(inter.shape[2]):
            inter[n,:,i] = np.linspace(inter[n,0,i],inter[n,-1,i],(X.shape[1]-1)*inflate+1)
    # return interpolated array
    return inter
