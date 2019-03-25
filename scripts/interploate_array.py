def interpolate(X,inflate):
    shape = (np.array(normal.shape) - 1 ) * inflate +1
    inter = np.zeros(shape=shape)
    inter[::inflate,::inflate,::inflate].shape

    inter[tuple(slice(None,None,inflate) for i in range(X.ndim))] = X
    
    for n_orig in range(X.shape[0]):
        n = n_orig * inflate
        for k_orig in range(X.shape[1]):
            k = k_orig*inflate
            for i_orig in range(X.shape[2]-1):
                i_new = i_orig*inflate
                j_new = (i_orig+1)*inflate
                inter[n,k][i_new:j_new+1] = np.linspace(inter[n,k][i_new],inter[n,k][j_new],inflate+1)

        for i in range(inter.shape[2]):
            inter[n,:,i] = np.linspace(inter[n,0,i],inter[n,-1,i],(X.shape[1]-1)*inflate+1)
        
    return inter
