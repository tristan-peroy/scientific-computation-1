"""
Code for Scientific Computation Project 1
Please add college id here
CID: 01854740
"""


#===== Code for Part 1=====#
from part1_utilities import * #do not modify

def method1(L_all,x):
    """
    First method for finding location of target in list containing
    M length-N sorted lists
    Input: L_all: A list of M length-N lists. Each element of
    L_all is a list of integers sorted in non-decreasing order.
    Example input for M=3, N=2: L_all = [[1,3],[2,4],[6,7]]

    """
    M = len(L_all)
    for i in range(M):
        ind = bsearch1(L_all[i],x)
        if ind != -1000:
            return((i,ind))

    return (-1000,-1000)




def method2(L_all,x,L_new = []):
    """Second method for finding location of target in list containing
    M length-N sorted lists
    Input: L_all: A list of M length-N lists. Each element of
    L_all is a list of integers sorted in non-decreasing order.
    Example input for M=3, N=2: L_all = [[1,3],[2,4],[6,7]]
    """

    if len(L_new)==0:
        M = len(L_all)
        N = len(L_all[0])
        L_temp = []
        for i in range(M):
            L_temp.append([])
            for j in range(N):
                L_temp[i].append((L_all[i][j],(i,j)))

        def func1(L_temp):
            M = len(L_temp)
            if M==1:
                return L_temp[0]
            elif M==2:
                return merge(L_temp[0],L_temp[1])
            else:
                return merge(func1(L_temp[:M//2]),func1(L_temp[M//2:]))

        L_new = func1(L_temp)

    ind = bsearch2(L_new,x)
    if ind==-1000:
        return (-1000,-1000),L_new
    else:
        return L_new[ind][1],L_new


def time_test(inputs=None):
    """Examine dependence of walltimes on M, N, and P for method1 and method2
        You may modify the input/output as needed.
    """

    #Add code here for part 1, question 2
    import numpy as np
    import time
    import matplotlib.pyplot as plt
    #initialize M,N,P
    M_length = 2**4
    N_length = 2**4
    P_length = [2**i for i in range(1,11)]
    
    #list to record relative runtimes of method1 and method2
    speeds1 = []
    speeds2 = []
    #iterate over increasing P
    for P in P_length:
            #generate random lists according to M,N,P
            L_all = [[np.random.randint(1,1000) for i in range(N_length)] for j in range(M_length)]
            targets = [0 for i in range(P)]
            t1 = time.time()
            #average over multiple runs
            for i in range(100):
                #for all P targets
                for target in targets:
                    y = method1(L_all,target)
            t2 = time.time()
            speed1 = (t2-t1)/100
            
            t3 = time.time()
            for i in range(100):
                L_new = []
                for target in targets:
                    y, L_new = method2(L_all,target,L_new)
            t4 = time.time()
            speed2 = (t4-t3)/100
            
            speeds1.append(speed1)
            speeds2.append(speed2)
    #plot both methods on the same graph
    plt.loglog(P_length, speeds1, base=2)
    plt.loglog(P_length, speeds2, base=2)
    plt.title("M=16,N=16, P increasing")
    plt.legend(["method1","method2"])
    plt.xlabel("log scale of P")
    plt.ylabel("log scale of worst-case runtime")
    plt.show()
    
    M_length = 1
    N_length = 3
    P_length = [2**i for i in range(1,11)]
    
    speeds1 = []
    speeds2 = []
    for P in P_length:
            L_all = [[np.random.randint(1,1000) for i in range(N_length)] for j in range(M_length)]
            targets = [0 for i in range(P)]
            t1 = time.time()
            for i in range(3000):
                for target in targets:
                    y = method1(L_all,target)
            t2 = time.time()
            speed1 = (t2-t1)/3000
            
            t3 = time.time()
            for i in range(3000):
                L_new = []
                for target in targets:
                    y, L_new = method2(L_all,target,L_new)
            t4 = time.time()
            speed2 = (t4-t3)/3000
            
            speeds1.append(speed1)
            speeds2.append(speed2)
    plt.loglog(P_length, speeds1, base=2)
    plt.loglog(P_length, speeds2, base=2)
    plt.title("M=1,N=3, P increasing")
    plt.legend(["method1","method2"])
    plt.xlabel("log scale of P")
    plt.ylabel("log scale of worst-case runtime")
    plt.show()



#===== Code for Part 2=====#

def findGene(L_in,L_p):
    """Find locations within adjacent strings (contained in input list,L_in)
    that contain patterns in input list L_p
    Input:
    L_in: A list containing two length-n strings
    L_p: A list containing p length-m strings

    Output:
    L_out: A length-p list whose ith element is a list of locations where the
    ith pattern has been found (see project description for further details)
    """
    #Size parameters
    n = len(L_in[0]) #length of a sequence
    p = len(L_p) #number of patterns
    m = len(L_p[0]) #length of pattern
    L_out = [[] for i in range(p)]

    #Add code here for part 2, question 1
    #using the Rabin Karp algorithm
    def heval(L,Base,Prime):
        """Convert list L to base-Base number mod Prime
        where Base specifies the base of L
        """
        f = 0
        for l in L[:-1]:
            f = Base*(l+f)
        h = (f + (L[-1])) % Prime
        return h
    def char2base4(S):
        """Convert gene test_sequence
        string to list of ints
        """
        c2b = {}
        c2b['A']=0
        c2b['C']=1
        c2b['G']=2
        c2b['T']=3
        L=[]
        for s in S:
            L.append(c2b[s])
        return L
    prime=997
    base = 4
    sequence0 = char2base4(L_in[0])
    sequence1 = char2base4(L_in[1])
    #convert both L_in sequences to base-4
    ind=0
    baseps = [char2base4(p) for p in L_p]
    hps = [heval(baseps[i],base,prime) for i in range(p)]
    hi1 = heval(sequence0[:m],base,prime)
    hi2 = heval(sequence1[:m],base,prime)
    #compute hash values for first m characters of both sequences
    for p, hp in enumerate(hps): #iterate over patterns
        if hi1==hp and hi2==hp: #Check if hash values match
            if sequence0[:m]==baseps[p] and sequence1[:m]==baseps[p]: #Character-by-character comparison
                L_out[p].append(ind) #Add index to output list

    bm = (4**m) % prime
    for ind in range(1,n-m+1): #iterate over remaining characters
    #Update rolling hash
        hi1 = (4*hi1 - int(sequence0[ind-1])*bm + int(sequence0[ind-1+m])) % prime
        hi2 = (4*hi2 - int(sequence1[ind-1])*bm + int(sequence1[ind-1+m])) % prime
        for p, hp in enumerate(hps):
            if hi1==hp and hi2==hp: #Check if hash values match
                if sequence0[ind:ind+m]==baseps[p] and sequence1[ind:ind+m]==baseps[p]: #Character-by-character comparison
                    L_out[p].append(ind)  #add index to output list


    return L_out


if __name__=='__main__':
    #Small example for part 2
    S1 = 'ATCGTACTAGTTATC'
    S2 = 'ATCTTAGTAGTCGTC'
    L_in = [S1,S2]
    L_p = ['ATC','AGT']
    out = findGene(L_in,L_p)

    #Large gene sequences
    infile1,infile2 = open("S1example.txt"), open("S2example.txt")
    S1,S2 = infile1.read(), infile2.read()
    infile1.close()
    infile2.close()
