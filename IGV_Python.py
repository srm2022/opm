# This is a generic Python implementation of IGV, proposed in the paper. Please customise it based on your requirements, e.g. by 
# changing the path to the data folder, the termination condition and the reported performance indicators.
# We customised its termination condition and reports differerntly for the purpose of comparing IGV with TS and IG on small and
# large datasets, as detailed in the paper. 
# Comments and questions may be emailed to seyed_r_mousavi@yahoo.com and new versions (if any) will be made available at
# https://github.com/srm2022/opm, as long as these accounts remain active. Thank you.

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
    print("main started")
    path = "drive/My Drive/OpM_LIB_2016/pmeds/"; # NOTE: Please change this path so it points to the folder (directory) on your machine where the data files are located.
    F = []  # this is the list fi in the paper (for client's facility)
    C = [[]] # this is the list of cj in the paper (for facility's clients)
    Cinx = [] # to indicate the index of each client in C
    nC = []
    
    P = [] # lists of centres (facilites)
    U = [] # list of non-centre vertices
    bestP = []
    runs = 10 # number of runs per instance
    
    file = open('IGV_Large.csv', 'a')
    # create the csv writer
    writer = csv.writer(file)
    data_list = "bestf_list", "duration_time_list"
    writer.writerow(data_list)
    
    for t in range(20, 40):
        print("t = ", t)
        filename = path + "pmed" + str(t + 1) + ".txt" # e.g. pmed21.txt
        prepareDandNandP(filename)
        for tmpk in range(0, 2):
            print("tmpk = ", tmpk)
            #overwrite p with n/4 or n/3:
            if (tmpk == 0):
                IGV_Python.p = int(IGV_Python.n/4) 
            else:
                IGV_Python.p = int(IGV_Python.n/3)
            print("actual p = ", IGV_Python.p)

            bestf_list = []
            duration_time_list = []

            for run in range(0,runs):
                print("run = ", run)
                start_time = time.time()
                IGV_Python.q = IGV_Python.n - IGV_Python.p
                P = []
                bestP = []
                U = []
      
                dim1, dim2 = (IGV_Python.n, IGV_Python.n) 
                C = [[0 for i in range(dim1)] for j in range(dim2)]
                IGV_Python.N = [[0 for i in range(dim1)] for j in range(dim2)]
                IGV_Python.Ninverse = [[0 for i in range(dim1)] for j in range(dim2)]
                IGV_Python.x = [False for i in range(dim1)]  
                dim1= IGV_Python.n
                Cinx = [0 for i in range(dim1)] 
                P = [0 for i in range(dim1)] 
                U= [0 for i in range(dim1)] 
                nC = [0 for i in range(dim1)] 
                F = [0 for i in range(dim1)]
                bestP = [0 for i in range(dim1)]
       
                fillin_N_and_Ninverse()
                g1 = g(1)
                f = initRandSol(P, U, IGV_Python.x, F, C, Cinx, nC)
                bestP = [P[i] for i in range(IGV_Python.p)]
                bestf = f
                radius = 1
                alpha = g1
                maxgap = 0
                duration = 10 # 10 seconds per each instance- please adjust the termination condition as you require
                
                while (time.time()-start_time < duration): # please adjust the termination condition based on your own requirements
                    # RLS1    
                    flgImproved = True
                    while (flgImproved and time.time()-start_time < duration): 
                        flgImproved = False

                        deltaf = 0
                        for d in range(0,radius):
                            delf= closefa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                            deltaf=deltaf+delf
                            IGV_Python.p=IGV_Python.p-1
                            IGV_Python.q=IGV_Python.q+1

                        for d in range(0,radius):
                            delf= openfa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                            deltaf=deltaf-delf
                            IGV_Python.p=IGV_Python.p+1
                            IGV_Python.q=IGV_Python.q-1

                        f += deltaf

                        if (deltaf > 0):
                            flgImproved = True
                            radius = 1
                            alpha = g1
                            if (f > bestf) :
                                for i in range(0,IGV_Python.p):
                                    bestP[i] = P[i]
                                maxgap = (maxgap * bestf + f - bestf)/f
                                bestf = f
                        pass

                    # RLS2
                    flgImproved = True
                    while (flgImproved and time.time()-start_time < duration): 
                        flgImproved = False

                        deltaf = 0
                        for d in range(0,radius):
                            delf = openfa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                            deltaf=deltaf-delf
                            IGV_Python.p=IGV_Python.p+1
                            IGV_Python.q=IGV_Python.q-1
                    
                        for d in range(0,radius):
                            delf = closefa(P, IGV_Python.p, U, IGV_Python.x, F, C, Cinx, nC, alpha)
                            deltaf=deltaf+delf
                            IGV_Python.p=IGV_Python.p-1
                            IGV_Python.q=IGV_Python.q+1                   

                        f += deltaf

                        if (deltaf > 0) :
                            flgImproved = True
                            radius = 1
                            alpha = g1
                            if (f > bestf) :
                                for i in range(0,IGV_Python.p):
                                    bestP[i] = P[i]
                                maxgap = (maxgap * bestf + f - bestf)/f
                                bestf = f

                    gap = (float(bestf) - f) / bestf
                    if (gap > maxgap):
                        maxgap = gap

                    if(maxgap == 0):
                        radius = 1
                    else:
                        radius += int ((1.0 * IGV_Python.p - 2) * gap/maxgap + 1)
                    if (radius > IGV_Python.p - 1):
                        radius = IGV_Python.p-1
                    alpha = g(radius)
                   
                equalsp=calculatef(bestP, IGV_Python.p)
                if (bestf !=equalsp): # to double-check the correctness of the algorithm
                        sys.exit("The final answer is incorrect!")
                bestf_list.append(bestf)
                duration_time_list.append(time.time()-start_time)
                           
            print("bestf_list = ", bestf_list, "duration_time_list = ", duration_time_list)
            data_list = bestf_list, duration_time_list
            writer.writerow(data_list)       
            print("---------------------------------------------------------")
           
    file.close()
    pass
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
    tmpabsdf=[]
    ret_val=0
    dim1= IGV_Python.n
    tmpabsdf = [0 for i in range(dim1)] 
    tmpq=IGV_Python.n-pp
    mindf = 99999999999 # use a sufficiently large number
    j_best_fa=-1
    for j in range(0,tmpq):
        ran_var=random.random()
        if(ran_var>alpha):
            continue
        fa=U[j]
        tmpabsdf[j]=0
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
                tmpabsdf[j_best_fa] =tmpabsdf[j_best_fa]+ IGV_Python.D[F[c]][c] - IGV_Python.D[fa][c]
                
        mindf = tmpabsdf[j_best_fa]
    #update
    ret_val = mindf
    best_fa = U[j_best_fa]
    P[pp] = best_fa
    pp=pp+1
    U[j_best_fa] = U[tmpq - 1]
    tmpq=tmpq-1
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
        
        pass
    return ret_val    
#----------------------
    
def closefa( P, pp, U, x,F, C, Cinx, nC, alpha):
    tmpabsdf=[]
    tmpq=IGV_Python.n-pp
    ret_val=0
    i_best_fa=-1
    maxdf=-1
    dim1= IGV_Python.n
    tmpabsdf = [0 for i in range(dim1)] 
    for i in range(0,pp):
        if random.random()>alpha:
            continue
        fa=P[i]
        tmpnc=nC[fa]
        tmpabsdf[i]=0
        for k in range(0,tmpnc):
            c=C[fa][k]
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
    tmpq=tmpq+1

    P[i_best_fa] = P[pp-1]
    pp=pp-1
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
    pass
#----------------------

def fillin_N_and_Ninverse():
    
    for c in range(0,IGV_Python.n):
        for j in range(0,IGV_Python.n):
            IGV_Python.N[c][j] = j
            #pass
        quicksort(c, 0, IGV_Python.n-1)
        pass
        
        #now Ninverse
    for c in range(0,IGV_Python.n):
        for j in range(0,IGV_Python.n):
            IGV_Python.Ninverse[c][IGV_Python.N[c][j]] = j
    pass
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
    pass
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
        if ( IGV_Python.x[fa]==False):
            U[k] = fa #or tmpU[j], the same
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
    pass
#----------------------

def nextpos(c, curpos):
    pos = 1 + curpos
    while(IGV_Python.x[IGV_Python.N[c][pos]]==False):
        pos=pos+1
    return pos
    pass
#----------------------

# this function is not needed for the large dataset but can be used to read instances of the small dataset
def read_n_and_D(theRawFilename):
    f = open(theRawFilename,'r')
    tmpS=f.read()
    posFrom = tmpS.index("n=") + 2
    posEnd = tmpS.index("m=") - 1
    IGV_Python.n = int(tmpS[posFrom: posEnd])
    posFrom = tmpS.index("m=") + 2
    posEnd = tmpS.index("p=") - 1
    tmpM = int((tmpS[posFrom: posEnd]))
    if (tmpM != IGV_Python.n):
        sys.exit("error: inconsistent n")
       
    posFrom = tmpS.index("p=") + 2
    posEnd = tmpS.index("client") - 1
    tmpP = int(tmpS[posFrom:posEnd])
    if (tmpP != IGV_Python.p):
        sys.exit("error: inconsistent p")
            
    dim1, dim2 = (IGV_Python.n, IGV_Python.n) 
    IGV_Python.D = [[-1 for i in range(dim1)] for j in range(dim2)] 

    posFrom = tmpS.index("table = ")
    posEnd = len(tmpS) - 1
    while (tmpS[posEnd] != '}') :
        posEnd=posEnd-1
    tmpS = tmpS[posFrom + 10: posEnd]

    for f in range(0,IGV_Python.n):
        posFrom = tmpS.index("{")
        posEnd = tmpS.index("}")
        tmpRow = tmpS[posFrom + 1: posEnd]
        tmpVals = tmpRow.split(",")
        if (f < IGV_Python.n - 1):
            tmpS = tmpS[posEnd + 2:]
        for c in range(0,IGV_Python.n):
            tmpfloat = float(tmpVals[c])
            IGV_Python.D[f][c] = int(tmpfloat)
            if (tmpfloat !=IGV_Python.D[f][c]):
                sys.exit("error: inconsistent weight")             
    pass
#----------------------

def prepareDandNandP(theRawFilename):
    tmpN = 0
    tmpM = 0
    tmpP = 0
    short_name=theRawFilename.split('\\')
    print(short_name)

    with open(theRawFilename) as f:
        lines = f.readlines()
    tmpParameters = lines[0].split(' ')
    tmpParameters=list(map(str.strip, tmpParameters))
    tmpParameters = list(filter(None,tmpParameters))  

    print(tmpParameters)
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
        #print(ijw)
        i=int(ijw[0])
        j=int(ijw[1])
        w=int(ijw[2])
        if ((i == j) or (w <= 0) or (i > tmpN) or (j > tmpN)):
            print("error: inconsistency in the input file!")
        else:
            IGV_Python.D[i - 1][j - 1] = IGV_Python.D[j - 1][i - 1] = w
    
    for k  in range(0,tmpN):
        for i  in range(0,tmpN):
            for j  in range(0,tmpN):
                if (IGV_Python.D[i][k] + IGV_Python.D[k][j] < IGV_Python.D[i][j]):
                    IGV_Python.D[i][j] = IGV_Python.D[i][k] + IGV_Python.D[k][j]
    
    IGV_Python.n = tmpN
    IGV_Python.p = tmpP

    pass
#----------------------

if __name__ == '__main__':
    main()
    pass
