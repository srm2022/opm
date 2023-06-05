
# This is a generic Python implementation of IGV. Please customise it as per 
# your own requirements, by changing its termination condition, reported 
# performance indicators, etc. We customised it differently for different 
# sections of the paper. Please also make sure the path points to your data 
# folder. To keep it simple, it only runs once, for 10 seconds, on each 
# instance of the large dataset and prints the best objective value and the 
# duration time (which is around 10 seconds, excluding the time spent on 
# reading input data). 
# Please note this code is meant to be simple and there is room for optimising it.
# Comments and questions are welcome to seyed_r_mousavi@yahoo.com and new 
# versions (if any) will be made available at # https://github.com/srm2022/opm, 
# as long as these accounts remain active. Thank you.

import sys
import math
import random
import time
from tokenize import Double
import numpy as np
import csv

class IGV_Python:
    n = 0
    p = 0
    q = 0
    D = [[]]
    N = [[]]
    Ninverse = [[]]
    x = [] # list of indicator variables xj. xi is true iff facility j is in P
    gamma = 0.6
    tau = 6

def main():
    path = "drive/My Drive/OpM_LIB_2016/pmeds/"; # NOTE: Please change this path so it points to the folder (directory) on your machine where the data files are located.
    F = []  # this is the list fi in the paper (for client's facility)
    C = [[]] # this is the list of cj in the paper (for facility's clients)
    Cinx = [] # to indicate the index of each client in C
    nC = []  
    g1 = g(1)
    for t in range(20, 40):
        filename = path + "pmed" + str(t + 1) + ".txt" # e.g. pmed21.txt
        prepareDandNandP(filename)
        
        for tmpk in range(0, 2):
            #overwrite p with n/4 or n/3:
            if (tmpk == 0):
                IGV_Python.p = int(IGV_Python.n/4) 
            else:
                IGV_Python.p = int(IGV_Python.n/3)
            print("pmed" + str(t + 1) + ", actual p = " + str(IGV_Python.p) + "\t", end='', flush=True)

            start_time = time.time()
            IGV_Python.q = IGV_Python.n - IGV_Python.p
      
            dim1 = IGV_Python.n
            P = [0 for i in range(dim1)] 
            bestP = [0 for i in range(dim1)]
            U= [0 for i in range(dim1)] 
            
            IGV_Python.N = [[0 for i in range(dim1)] for j in range(dim1)]
            IGV_Python.Ninverse = [[0 for i in range(dim1)] for j in range(dim1)]
            fillin_N_and_Ninverse()
            
            IGV_Python.x = [False for i in range(dim1)]  
            F = [0 for i in range(dim1)]
            C = [[0 for i in range(dim1)] for j in range(dim1)]
            Cinx = [0 for i in range(dim1)] 
            nC = [0 for i in range(dim1)] 

            f = initRandSol(P, U, IGV_Python.x, F, C, Cinx, nC)

            bestP = [P[i] for i in range(IGV_Python.p)]
            worstf = bestf = f
            radius = 1
            alpha = g1
            duration = 10 # 10 seconds per instance
                
            while (time.time()-start_time < duration): # please adjust the termination condition based on your own requirements
                # RLS1    
                flgImproved = True
                while (flgImproved and time.time()-start_time < duration): 
                    flgImproved = False

                    deltaf = 0
                    for d in range(0,radius):
                        delf = closefa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                        deltaf = deltaf + delf
                        IGV_Python.p = IGV_Python.p - 1
                        IGV_Python.q = IGV_Python.q + 1

                    for d in range(0,radius):
                        delf = openfa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                        deltaf = deltaf - delf
                        IGV_Python.p = IGV_Python.p + 1
                        IGV_Python.q = IGV_Python.q - 1

                    f += deltaf

                    if (deltaf > 0):
                        flgImproved = True
                        radius = 1
                        alpha = g1
                        if (f > bestf) :
                            for i in range(0, IGV_Python.p):
                                bestP[i] = P[i]
                            bestf = f
                    if( f < worstf):
                      worstf = f

                # RLS2
                flgImproved = True
                while (flgImproved and time.time()-start_time < duration): 
                    flgImproved = False

                    deltaf = 0
                    for d in range(0, radius):
                        delf = openfa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                        deltaf = deltaf - delf
                        IGV_Python.p = IGV_Python.p + 1
                        IGV_Python.q = IGV_Python.q - 1
                   
                    for d in range(0, radius):
                        delf = closefa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                        deltaf = deltaf + delf
                        IGV_Python.p = IGV_Python.p - 1
                        IGV_Python.q = IGV_Python.q + 1                   

                    f += deltaf

                    if (deltaf > 0) :
                        flgImproved = True
                        radius = 1
                        alpha = g1
                        if (f > bestf) :
                            for i in range(0, IGV_Python.p):
                                bestP[i] = P[i]
                            bestf = f
                    if( f < worstf):
                      worstf = f

                if (worstf == bestf):
                    gap = 0
                else:
                    gap = (float(bestf) - f) / (bestf - worstf)

                radius += int ((1.0 * IGV_Python.p - 2) * gap + 1)
                if (radius > IGV_Python.p - 1):
                    radius = IGV_Python.p - 1

                alpha = g(radius)
                   
            equalsp=calculatef(bestP, IGV_Python.p)
            if (bestf !=equalsp): # to double-check the correctness of the algorithm
                sys.exit("The final answer is incorrect!")
            print(bestf, time.time()-start_time, flush=True)                           
#----------------------

def calculatef(pP,pp):
    ret_val=0
    for c in range(0,IGV_Python.n):
        tmpmin = math.inf
        for i in range(0,pp):
            if(IGV_Python.D[pP[i]][c] < tmpmin):
                tmpmin = IGV_Python.D[pP[i]][c]
        ret_val+=tmpmin
    
    return ret_val
#----------------------

def openfa( P, pp, U, x,F, C, Cinx, nC, alpha):
    ret_val = 0
    dim1= IGV_Python.n
    tmpabsdf = [0 for i in range(dim1)] 
    tmpq=IGV_Python.n-pp
    mindf = 999999999 # use a sufficiently large number
    j_best_fa = -1
    for j in range(0,tmpq):
        ran_var = random.random()
        if(ran_var > alpha):
            continue
        fa = U[j]
        tmpabsdf[j] = 0
        for c in range(0,IGV_Python.n):
            if (IGV_Python.D[fa][c] < IGV_Python.D[F[c]][c]):
                tmpabsdf[j] = tmpabsdf[j]+IGV_Python.D[F[c]][c] - IGV_Python.D[fa][c]
        if (tmpabsdf[j] < mindf):
            mindf = tmpabsdf[j]
            j_best_fa = j  
    if (j_best_fa == -1):
        j_best_fa = int (random.random() * tmpq)

        fa = U[j_best_fa]
        tmpabsdf[j_best_fa] = 0
        for c in range(0,IGV_Python.n):
            #calculate contributation of c to tmpdf
            if (IGV_Python.D[fa][c] < IGV_Python.D[F[c]][c]):
                tmpabsdf[j_best_fa] += IGV_Python.D[F[c]][c] - IGV_Python.D[fa][c]
                
        mindf = tmpabsdf[j_best_fa]
    #update
    ret_val = mindf
    best_fa = U[j_best_fa]
    P[pp] = best_fa
    pp += 1
    U[j_best_fa] = U[tmpq - 1]
    tmpq -= 1
    IGV_Python.x[best_fa] = True

    nC[best_fa] = 0  # //CRUCIAL     
    for c in range(0,IGV_Python.n):
        if (IGV_Python.Ninverse[c][best_fa] < IGV_Python.Ninverse[c][F[c]]):  # Please note that we CANNOT use if (D[best_fa][c] < D[F[c]][c])
            #remove from old list
            oldfa = F[c]
            oldinx = Cinx[c]
            tmpVlast = C[oldfa][nC[oldfa] - 1]
            C[oldfa][oldinx] = tmpVlast
            Cinx[tmpVlast] = oldinx
            nC[oldfa]=nC[oldfa]-1

            #add to new list
            F[c] = best_fa
            newinx = nC[best_fa]
            C[best_fa][newinx] = c
            Cinx[c] = newinx
            nC[best_fa]=nC[best_fa]+1
        
    return ret_val    
#----------------------
    
def closefa( P, pp, U, x,F, C, Cinx, nC, alpha):
    tmpq = IGV_Python.n-pp
    ret_val = 0
    i_best_fa = -1
    maxdf = -1
    dim1= IGV_Python.n
    tmpabsdf = [0 for i in range(dim1)] 
    for i in range(0,pp):
        if random.random()>alpha:
            continue
        fa = P[i]
        tmpnc = nC[fa]
        tmpabsdf[i] = 0
        for k in range(0,tmpnc):
            c = C[fa][k]
            tmpnextpos = nextpos(c, IGV_Python.Ninverse[c][F[c]])
            nextminfa = IGV_Python.N[c][tmpnextpos]               
            mind = IGV_Python.D[nextminfa][c]
                
            tmpabsdf[i] += mind - IGV_Python.D[fa][c]
        if (tmpabsdf[i] > maxdf):
            maxdf = tmpabsdf[i]
            i_best_fa = i

    if (i_best_fa == -1):
        i_best_fa = int(random.random() * pp)
        fa = P[i_best_fa]
        tmpnc = nC[fa]
        #calculate the amount of change in f if fa dropped
        tmpabsdf[i_best_fa] = 0
        for k in range(0,tmpnc):
            c = C[fa][k]
            tmpnextpos = nextpos(c, IGV_Python.Ninverse[c][F[c]])
            nextminfa = IGV_Python.N[c][tmpnextpos];                
            mind = IGV_Python.D[nextminfa][c]

            tmpabsdf[i_best_fa] += mind - IGV_Python.D[fa][c]     
        maxdf = tmpabsdf[i_best_fa]
    #update info
    ret_val = maxdf

    best_fa = P[i_best_fa]
    U[tmpq] = best_fa
    tmpq += 1

    P[i_best_fa] = P[pp-1]
    pp -= 1
    IGV_Python.x[best_fa] = False

    #update F and C as well
    tmpnc = nC[best_fa]
    for k in range(0,tmpnc): 
        c = C[best_fa][k]
        tmpnextpos = nextpos(c, IGV_Python.Ninverse[c][F[c]])
        tmpminfa = IGV_Python.N[c][tmpnextpos]
        F[c] = tmpminfa
        tmpinx = nC[tmpminfa]
        C[tmpminfa][tmpinx] = c
        Cinx[c] = tmpinx
        nC[tmpminfa]=nC[tmpminfa]+1
    nC[best_fa] = 0 # //not required, just for more reliability
    return ret_val         
#----------------------

def g(radius):
    return IGV_Python.gamma *pow(2, IGV_Python.tau*(1-radius)) + (1-IGV_Python.gamma) * random.random()
#----------------------

def fillin_N_and_Ninverse():
    
    for c in range(0,IGV_Python.n):
        for j in range(0,IGV_Python.n):
            IGV_Python.N[c][j] = j
        quicksort(c, 0, IGV_Python.n-1)
        
        #now Ninverse
    for c in range(0,IGV_Python.n):
        for j in range(0,IGV_Python.n):
            IGV_Python.Ninverse[c][IGV_Python.N[c][j]] = j
#----------------------

def quicksort(c,From,to):    
    if(From < to):
        pivotitem = IGV_Python.D[IGV_Python.N[c][From]][c]
        j = From
        for i in range(From+1,to+1):
            if(IGV_Python.D[IGV_Python.N[c][i]][c] < pivotitem):
                j=j+1
                tmp = IGV_Python.N[c][i]
                IGV_Python.N[c][i] = IGV_Python.N[c][j]
                IGV_Python.N[c][j] = tmp
                
        tmp = IGV_Python.N[c][j]
        IGV_Python.N[c][j] = IGV_Python.N[c][From]
        IGV_Python.N[c][From] = tmp
        quicksort(c, From, j - 1)
        quicksort(c, j + 1, to)
#----------------------

def initRandSol(P, U, x, F, C, Cinx, nC):
    tmpPbyN = IGV_Python.p * 1.0 / IGV_Python.n
    k = 0
    fa = int (random.random() * IGV_Python.n) # //random starting point
    while (k < IGV_Python.p):
        if (IGV_Python.x[fa]==False):
            if (random.random() < tmpPbyN):
                P[k] = fa
                IGV_Python.x[fa] = True
                k=k+1
            
        fa = (fa + 1) % IGV_Python.n
    k = 0

    for fa in range(0,IGV_Python.n):
        if ( IGV_Python.x[fa] == False):
            U[k] = fa 
            k=k+1
        
    tmpf = 0
    for c in range(0,IGV_Python.n):
        minpos = nextpos(c, -1)
      
        minfa = IGV_Python.N[c][minpos]          
        mind = IGV_Python.D[minfa][c]

        F[c] = minfa
        tmpinx = nC[minfa]
        C[minfa][tmpinx] = c
        Cinx[c] = tmpinx
        nC[minfa]=nC[minfa]+1
        tmpf += mind
    return tmpf
#----------------------

def nextpos(c, curpos):
    pos = 1 + curpos
    while(IGV_Python.x[IGV_Python.N[c][pos]]==False):
        pos=pos+1
    return pos
#----------------------

def prepareDandNandP(theRawFilename):
    tmpN = 0
    tmpM = 0
    tmpP = 0
    short_name=theRawFilename.split('\\')

    with open(theRawFilename) as f:
        lines = f.readlines()
    tmpParameters = lines[0].split(' ')
    tmpParameters=list(map(str.strip, tmpParameters))
    tmpParameters = list(filter(None,tmpParameters))  

    tmpN = int(tmpParameters[0])
    tmpM = (int)(tmpParameters[1])
    tmpP = (int)(tmpParameters[2])

    IGV_Python.D=[[]]
    IGV_Python.D = [[0 for i in range(tmpN)] for j in range(tmpN)]
    for i in range(0,tmpN):
        IGV_Python.D[i][i]=0
        for j in range(i+1,tmpN):
            IGV_Python.D[i][j] = IGV_Python.D[j][i] = math.inf
        
    for i in range(1,len(lines)):
        ijw = lines[i].split(' ')
        ijw=list(map(str.strip, ijw))
        ijw = list(filter(None,ijw))
        i=int(ijw[0])
        j=int(ijw[1])
        w=int(ijw[2])
        if ((i == j) or (w <= 0) or (i > tmpN) or (j > tmpN)):
            print("error: inconsistency in the input file!")
        else:
            IGV_Python.D[i - 1][j - 1] = IGV_Python.D[j - 1][i - 1] = w
    
    #make the complete graph by finding the shortest paths via Floyd–Warshall alg.
    for k  in range(0,tmpN):
        for i  in range(0,tmpN):
            for j  in range(0,tmpN):
                if (IGV_Python.D[i][k] + IGV_Python.D[k][j] < IGV_Python.D[i][j]):
                    IGV_Python.D[i][j] = IGV_Python.D[i][k] + IGV_Python.D[k][j]
    
    IGV_Python.n = tmpN
    IGV_Python.p = tmpP

#----------------------

if __name__ == '__main__':
    main()
