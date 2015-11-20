from scipy import stats
import scipy
import numpy as np

'''
Module for studying fixation and extinction properties of the neutral growth model.
'''

def firstbottleneck(N=10, Ndiv=5, initialfrequency=0.05):
    '''
    Set up the probabilities p for mutation frequencies during the first bottleneck after a mutation
    arises. The frequency before the bottleneck will generally be low since the mutation arises in one of several
    individuals, while after the bottleneck it will be in the set {0,1/N,2/N,...,N-1/N,1}
    p[i] contains the The probabilities of frequency i/N is given by the hypergeometric distribution f(i;N*2**Ndiv,initialfrequency*2**Ndiv,N)
    
    '''
    p=np.zeros(N+1)
    
    for i in range(N+1):
        p[i]=stats.hypergeom.pmf(i,N*2**Ndiv,initialfrequency*N*2**Ndiv,N)
        #print i, N*2**Ndiv, initialfrequency*N*2**Ndiv

    return p


def bottleneckmatrix(N=10, Ndiv=5):
    '''
    Set up the stochastic matrix X modelling the Markov process occuring during
    a bottleneck event where N out of 2^Ndiv*N individuals are sampled. It is assumed 
    that the mutation was already there during the preceeding bottleneck and so it is limited to
    N+1 different frequencies: {0,1/N,2/N,...,N-1/N,1}

    Element X[i,j] describes of a mutation with frequency (j-1)/N before the bottleneck 
    having frequency (i-1)/N after the bottleneck and this probability is given by the
    hypergeometric distribution (http://en.wikipedia.org/wiki/Hypergeometric_distribution)
    with parameters f(i;N*2**Ndiv,j*2**Ndiv,N)
    '''
    X=np.zeros((N+1,N+1))*np.nan

    for i in range(N+1):
        for j in range(N+1):
            X[i,j]=stats.hypergeom.pmf(i,N*2**Ndiv,j*2**Ndiv,N)
    
    #special case where frequency is 1
    X[:,0]=np.zeros(N+1)
    X[0,0]=1.0

    return X

def expected_time_to_extinction_and_fixation(N=10,Ndiv=5,Nbottle=1000):
    '''
    Calculate expected time to extinction and fixation for a weighted average
    of mutations using default parameters.
    '''
    
    #set up weighted average of frequencies after first bottleneck
    # mutation happening during n-th cell division appear in frequency 1/2^n*N
    x01 = 10*firstbottleneck(initialfrequency=1.0/20.0)
    x02 = 20*firstbottleneck(initialfrequency=1.0/40.0)
    x03 = 40*firstbottleneck(initialfrequency=1.0/80.0)
    x04 = 80*firstbottleneck(initialfrequency=1.0/160.0)
    x05 = 160*firstbottleneck(initialfrequency=1.0/320.0)

    #weighted average, give expected distribution of mutation frequencies after first bottleneck
    #after mutation occured
    x0 = np.atleast_2d((x01+x02+x03+x04+x05)/310)
    x0.shape = (11,1)
    #segsum=sum(x0[range(1,10),:])[-1] #expected frequency of mutations still segregating

    # Markov matrix representing the effect of one bottleneck event
    X = bottleneckmatrix()

    #simulate bottlenecks
    for i in range(Nbottle-1):
        x0=np.column_stack((x0,scipy.dot(X,x0[:,-1])))
        segsum=sum(x0[range(1,10),:])[-1] #expected frequency of mutations still segregating

    #add zeros as first column for scipy.diff below to count the first bottleneck 
    x0 = np.column_stack((np.zeros((11,1)),x0))
    
    #print results
    print "Ndiv:                         ",Ndiv
    print "N:                    ",N
    print "Nbottle:              ",Nbottle
    print "Rate of extinction:      ",x0[0,-1]
    print "Mean time to extinction: ",scipy.sum(scipy.diff(x0[0,:])*np.array(range(1,x0.shape[1])))/x0[0,-1]
    print "Rate of fixation:        ",x0[-1,-1]
    print "Mean time to fixation:   ",scipy.sum(scipy.diff(x0[-1,:])*np.array(range(1,x0.shape[1])))/x0[-1,-1]
    print "Still segregating:       ",segsum

expected_time_to_extinction_and_fixation()