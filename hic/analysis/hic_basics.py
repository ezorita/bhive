from scipy import sparse
import scipy
import scipy.linalg
import numpy as np
import cooler
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


  
def scale_resolution(matrix, factor):
    new_row = matrix.row // factor
    new_col = matrix.col // factor
    bin = (matrix.shape[0] // factor) * new_col + new_row
    unique, bin = np.unique(bin, return_inverse=True)
    sum = np.bincount(bin, weights=matrix.data)
    new_col = unique // (matrix.shape[0] // factor)
    new_row = unique - new_col * (matrix.shape[0] // factor)

    return scipy.sparse.coo_matrix((sum, (new_row, new_col)))
        

def cpb_step(m):
    # Matrix cpb balancing.
    its = m.sum(axis=1)        # Total interaction signal
    itc = (m>0).sum(axis=1) # Number of interactions
    cpb = its/(itc.astype('f'))
    cpb = np.asarray(cpb).reshape(cpb.shape[0],)
    # Normalize CPB
    cpb_notnan = cpb[np.where(np.logical_not(np.isnan(cpb)))]
    cpb = cpb/np.max(cpb_notnan)
    cpb[np.where(np.isnan(cpb))] = 1.0
    # Make diagonal matrix containing cpb values
    d = sparse.coo_matrix(m.shape, dtype=np.float)
    d.setdiag(np.sqrt(1.0/cpb))
    m = d*m*d
    # Recompute variance
    #    cpb = np.asarray(m.sum(axis=1)/itc.astype('f')).reshape(itc.shape[0],)
    #    cpb_notnan = cpb[np.where(np.logical_not(np.isnan(cpb)))]
    return m #, np.var(cpb_notnan/np.mean(cpb_notnan))

def cpb_var(m):
    # Matrix cpb balancing.
    its = m.sum(axis=1)        # Total interaction signal
    itc = (m>0).sum(axis=1) # Number of interactions
    cpb = its/(itc.astype('f'))
    cpb = np.asarray(cpb).reshape(cpb.shape[0],)
    # Normalize CPB
    cpb_notnan = cpb[np.where(np.logical_not(np.isnan(cpb)))]
    cpb = cpb/np.max(cpb_notnan)
    return np.var(cpb_notnan/np.mean(cpb_notnan))

def cpb(matrix, n_iter):
    '''
    Contact Probability Balancing algorithm.
    Returns CPB-balanced version of matrix after applying n_iter iterations.
    '''
    m = sparse.coo_matrix(matrix, dtype=float)
    #cvar = np.array([0.0]*(n_iter+1))
    #cvar[0] = cpb_var(m)
    for i in xrange(1,n_iter+1):
        m = cpb_step(m)
        
    return m #, cvar

def oe_apply(m, coeffs):
    '''
    Applies O/E normalization given a sparse/coo_matrix and previously computed coefficients.
    '''
    matrix = sparse.coo_matrix(m, dtype=float)
    for i in range(0, len(matrix.data)):
        matrix.data[i] = matrix.data[i] / coeffs[abs(matrix.col[i]-matrix.row[i])]

    return matrix
    
def oe_coeffs_chrom(cfile, ranges):
    '''
    Returns a dictionary with O/E normalization coefficients for each chromosome.
    Dictionary keys are chromosome names.
    '''
    if type(ranges) is not list:
        ranges = [ranges,]
    pr_int = {}
    # Compute observed interactions for selected chromosomes
    chrlist = []
    chrbins = []
    for cname in ranges:
        chrname = cname.split(':')[0]
        chrlist.append(chrname)
        chrbins.append(cfile.extent(str(chrname)))
 
    for c,bins in zip(chrlist,chrbins):
        max_int_dist = bins[1] - bins[0]
        # ObI: Observed Interactions
        # PrI: Probability of Interaction
        obs_int = [0]*(max_int_dist+1)
        pr_int[c] = [0.0]*(max_int_dist+1)

        # Compute observed interactions for each chromosomes
        m = cfile.matrix(balance=False, sparse=True).fetch(c)
        for i,j,v in zip(m.row, m.col, m.data):
            obs_int[abs(i-j)] += v

        for int_dist in range(0,max_int_dist+1):
            if not obs_int[int_dist]: continue
            # PsI: Possible interactions
            possible_int = max_int_dist-int_dist
            # Probability:
            pr_int[c][int_dist] = obs_int[int_dist] / float(possible_int)

    return pr_int

def oe_coeffs_gw(cfile):
    '''
    Returns an array containing the O/E normalization coefficient for all
    interaction distances.
    '''
    # Compute observed interactions for selected chromosomes
    chrlist = []
    chrbins = []
    for cname in cfile.chroms()[:]['name']:
        chrname = cname.split(':')[0]
        chrlist.append(chrname)
        chrbins.append(cfile.extent(str(chrname)))
            
    # Max interaction distance
    chrsz = sorted([x[1]-x[0] for x in chrbins], reverse = 1)
    max_int_dist = chrsz[0]-1
    
    # ObI: Observed Interactions
    # PrI: Probability of Interaction
    obs_int = [0]*(max_int_dist+1)
    pr_int = [0.0]*(max_int_dist+1)

    # Compute observed interactions for all chromosomes
    for c in chrlist:
        m = cfile.matrix(balance=False, sparse=True).fetch(c)
        for i,j,v in zip(m.row, m.col, m.data):
            obs_int[abs(i-j)] += v

    # Compute probability of interaction
    nchr = len(chrsz)
    cumsz = sum(chrsz)

    for int_dist in xrange(0,max_int_dist+1):
        if obs_int[int_dist] == 0: continue
        # Update chromosome count
        while chrsz[nchr-1] <= int_dist:
            nchr -= 1
            cumsz -= chrsz[nchr]
        # PsI: Possible interactions
        possible_int = cumsz-nchr*int_dist
        # Probability:
        pr_int[int_dist] = obs_int[int_dist] / float(possible_int)

    pr_int = np.array(pr_int)
    pr_int[np.where(pr_int == 0)] = np.min(pr_int[np.where(pr_int > 0)])
            
    return pr_int

def eigs_ab_score(cor_mat, ab_percentile=20, cormat_percentile=20):
    '''
    Inputs:
    sparse.coo_matrix cor_mat : correlation matrix
    ab_percentile: percentile of reference A/B compartments
    Output:
    returns AB score for all bins in cor_mat
    '''
    # Block eigenvector computation (of correlation matrix):
    eigv = np.empty(cor_mat.shape[0], dtype='f')
    ab = np.empty(cor_mat.shape[0], dtype='f')
    ab_score = np.empty(cor_mat.shape[0], dtype='f')
    eigv[:] = np.log(-1)
    ab[:] = 0
    ab_score[:] = 0
        
    notnan = []
    for i in xrange(0,cor_mat.shape[0]):
        if not np.isnan(cor_mat[i,:]).all() and not (cor_mat[i,:] == 0).all():
            notnan.append(i)

    notnan = np.array(notnan)
    # Center matrix to median within cormat_percentiles
    cm = cor_mat[notnan,:][:,notnan]
    offdiag = np.extract(1-np.eye(cm.shape[0]),cm)
    # Limit value distribution to cormat_percentiles
    # Update: Compute percentiles of negative and positive values independently
    low = np.percentile(offdiag[np.where(offdiag < 0)],cormat_percentile)
    up  = np.percentile(offdiag[np.where(offdiag > 0)],100-cormat_percentile)
    cor_mat[np.where(cor_mat<low)] = low
    cor_mat[np.where(cor_mat>up )] = up

    # Center to median (maybe not a good idea...)
    med = np.median(offdiag)
    cor_mat = cor_mat - med
    low -= med
    up -= med
    
    # Scale to -1, 1
    cor_mat[np.where(cor_mat<0)] = cor_mat[np.where(cor_mat<0)]/np.abs(low)
    cor_mat[np.where(cor_mat>0)] = cor_mat[np.where(cor_mat>0)]/np.abs(up)

    # Erase diagonal
    np.fill_diagonal(cor_mat,0)
    cor_mat[np.where(np.isnan(cor_mat))] = 0
    
    # Compute eigenvectors
    N = cor_mat.shape[0]
    w, v = scipy.linalg.eigh(cor_mat, eigvals=(N-1,N-1))
    v = v[:,np.argmax(w)]
    v = v/np.max(np.abs(v))

    # Compute A/B scores
    a_set = np.where(v > np.percentile(v,100-ab_percentile))
    b_set = np.where(v < np.percentile(v,ab_percentile))
    for i in xrange(0,cor_mat.shape[0]):
        a_score = np.sum(cor_mat[i,:][a_set])/float(len(a_set))
        b_score = np.sum(cor_mat[i,:][b_set])/float(len(b_set))
        a_mod   = np.sum(np.abs(cor_mat[i,:][a_set]))/float(len(a_set))
        b_mod   = np.sum(np.abs(cor_mat[i,:][b_set]))/float(len(b_set))
        ab_score[i] = (a_score-b_score)/(a_mod+b_mod)*100.0

    # Store data
    eigv[notnan] = v

    return eigv, ab_score

