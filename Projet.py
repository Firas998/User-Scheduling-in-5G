import copy
import numpy as np
from functools import cmp_to_key
import matplotlib.pyplot as plt
from scipy.optimize import linprog
import random
from statistics import mean
from copy import deepcopy
from time import perf_counter 



def parsefile(path):
    datar=[]
    datap=[]
    lines=open(path,'r').readlines()
    lines = [line.rstrip() for line in lines]
    #First we extract N,M,K,P
    N=int(float(lines[0]))
    M=int(float(lines[1]))
    K=int(float(lines[2]))
    P=int(float(lines[3]))
    indice=4
    #Second, we create an array of N arrays ( one for each channel, which contains all powers )
    for i in range(N):
        curdata=[]
        for e in range (K):
            curline=lines[indice].split()
            curline  = [ int(float(x)) for x in curline ]
            indice+=1
            curdata.extend(curline)
        datap.append(curdata)
    #Last, we create an array of N arrays ( one for each channel, which contains all rates)
    for i in range(N):
        curdata=[]
        for e in range (K):
            curline=lines[indice].split()
            curline  = [ int(float(x)) for x in curline ]
            indice+=1
            curdata.extend(curline)
        datar.append(curdata)
    return (datap,datar,P)
  
def preprocessq2(fulldata,powerbudget):
    n=len(fulldata)
    #We create a table that contains the minimum power of each channel
    P=[] 
    for i in range(n):
        values=[item[0] for item in fulldata[i]]
        P.append(min(values))
        
    somme=sum(P) #The sum of all minimum powers
    for i in range(n):
        cursomme=somme-P[i]
        for e in range(len(values)-1, -1, -1):
            if fulldata[i][e][0]+cursomme>powerbudget:
                fulldata[i].pop(e)

def compare(item1, item2): #This is the comparison function we will use in the sorting of all pairs for removing IP-dominated terms
    if item1[0]<item2[0]:
        return -1
    elif item1[0] > item2[0]:
        return 1
    elif item1[1]<item2[1]:
        return 1
    else:
      return -1
def preprocessq3(data):
    n=len(data)
    for i in range(n):
        values=data[i]
        if len(values)==0:
            continue
        values=sorted(values, key=cmp_to_key(compare))
        maximum=values[0][1]
        valuescopy=[[values[0][0],maximum]]
        for e in range(1,len(values)):
            if values[e][1]>maximum:
                maximum=values[e][1]
                valuescopy.append([values[e][0],maximum])
        data[i]=valuescopy

#We define a function that gives us the slope between two different pairs
slope = lambda item1,item2 : (item2[1]-item1[1])/float(item2[0]-item1[0])

def preprocessq4(data):
  for i in range(len(data)):
      Values=data[i]
      if len(Values)<2:
          continue
      Resultat=[Values[0],Values[1]]
      for z in range(2,len(Values)):
          Nov=Values[z]
          while len(Resultat)>=2 and slope(Resultat[-2],Resultat[-1])<=slope(Resultat[-2],Nov):
            Resultat.pop()
          Resultat.append(Nov)
      data[i]=Resultat
          
def greedy (fulldata,p,startingpoint):
    start= perf_counter()  
    fulldata=fulldata[startingpoint:]
    n=len(fulldata)
    remainpower=p
    result=0
    #First we process the first elements for each channel
    for i in range(n):
        if len(fulldata[i])==0:
            continue
        remainpower-=fulldata[i][0][0]
        result+=fulldata[i][0][1]
    #Now we restructure the pairs in every channel such that we have access to increment, increment efficiency, and increment power
    for i in range(n):
        values=fulldata[i]
        for e in range(1,len(values)):
            pair=values[e]
            prev=values[e-1]
            values[e]=[pair[0],pair[1],(pair[1]-prev[1])/(pair[0]-prev[0]),pair[1]-prev[1],pair[0]-prev[0]]
    #We put all pairs in a single table and sort them
    allvalues=[]
    for table in fulldata:
        allvalues.extend(table[1:])
    allvalues = sorted(allvalues, key=lambda x:-x[2])
    
    #Now we apply the greedy algorithm, taking one by one the pairs in descending order of incremental efficiency
    increment=0
    add=0
    testenough=False #This boolean determines if we left the loop because of an insufficient power or not
    while len(allvalues)>0:
        pair=allvalues.pop(0)
        increment=pair[4]
        add=pair[3]
        if remainpower<increment:
            testenough=True
            break
        remainpower-=increment
        result+=pair[3]
    intsol=result
    #In the case insufficient power was to be held accountable for the end of the loop we add the surplus
    if testenough:
        surplus=remainpower/increment
        result+=surplus*add
    end=perf_counter()  
    return [intsol,result,p-remainpower,p,end-start]


def LPsolveprocessed(fulldata,powerbudget):
    starttime=perf_counter()  
    N=len(fulldata)
    coeffs=[]
    totallen=0
    for i in range(N):
        Values=fulldata[i]
        totallen+=len(Values)
        coeffs.extend([-Values[e][1] for e in range(len(Values)) ])
    if totallen==0:
        return [0,0]
    #We build A_eq and B_eq which serve to ensure the constraint that the sum of the binary values is 1 for each channel
    A_eq=[]
    toadd=[]
    remain=totallen
    start=0
    for n in range(N):
        toadd=[0 for e in range(start)]
        Values=fulldata[n]
        remain-=len(Values)
        start+=len(Values)
        toadd.extend([1 for e in range(len(Values))])
        toadd.extend([0 for e in range(remain)])
        A_eq.append(toadd)
    B_eq=[1 for i in range (N)]
    
    #We build A_ub and B_ub which serve to ensure the constraint that we do not surpass the powerbudget
    powercheck=[]
    for i in range(N):
        Values=fulldata[i]
        powercheck.extend([Values[e][0] for e in range(len(Values))])
    A_ub=[powercheck]
    B_ub=[[powerbudget]]
    
    positif=[(0,1) for i in range(totallen)]
    lpsolved=linprog(c=coeffs,A_ub=A_ub,b_ub=B_ub,A_eq=A_eq,b_eq=B_eq,bounds=positif)
    endtime= perf_counter()  
    
    return [-lpsolved.fun,endtime-starttime]
        
def dynamic1(fulldata,powerbudget):
    start=perf_counter()
    N=len(fulldata)
    Resultat=[0 for col in range(powerbudget)]
    Values=fulldata[0]
    #First we check if a a channel has no potential pairs, in which case there is no solution
    for i in range(N):
        if len(fulldata[i])==0:
            return [0,0,0]
    #We fill the table for the first channel
    for i in range(powerbudget):
        liste=[pair[1] for pair in Values if pair[0]<=(i+1)]
        if len(liste)!=0:
            Resultat[i]=max(liste)
    #We then apply the dynamic programming formula
    for i in range(1,N):
        Values=fulldata[i]
        for power in  range(powerbudget-1, -1, -1): #We iterate backwards so as not to cause any conflict between different values
            maximum=0
            for pair in Values:
                if pair[0]<=(power) and Resultat[power-pair[0]]>0:
                    maximum=max(maximum,pair[1]+Resultat[power-pair[0]])
            Resultat[power]=maximum
    #We check for the powerbudget used, which is the index of the entry in the table with the same rate as the highest rate
    for i in range(len(Resultat)):
        if Resultat[i]==Resultat[-1]:
            end=perf_counter()
            return [Resultat[-1],i,end-start]

def dynamic2(fulldata,powerbudget):
    
    start=perf_counter()
    N=len(fulldata)
    #We first we check that all channels have at least one pair, otherwise there is no feasible solution
    for i in range(N):
        if len(fulldata[i])==0:
            return [0,0]
    U=0
    for i in range(len(fulldata)):
        liste=[item[1] for item in fulldata[i]]
        U+=max(liste)    
    U+=1
    Resultat=[0 for i in range(U)]
    #We start by filling the table for the first channel
    for i in range(U):
        liste=[item[0] for item in fulldata[0] if item[1]==i]
        if (len(liste)!=0):
            Resultat[i]=min(liste)
    #We apply the dynamic programming formula
    for i in range(1,N):
        Values=fulldata[i]
        for rate in range(U-1, -1, -1):  #We iterate backwards so as not to cause any conflict between values   
            liste=[]
            for pair in Values:
                if pair[1]<=(U+1) and Resultat[rate-pair[1]]!=0:
                    liste.append(Resultat[rate-pair[1]]+pair[0])
            if (len(liste)!=0):
                Resultat[rate]=min(liste)
    for rate in range(U-1, -1, -1): #We iterate backwards to find the highest possible rate, which is the first rate encountered backwards that fulfills our constraints
        if (Resultat[rate]<=powerbudget and Resultat[rate]!=0):
            end=perf_counter()
            return [rate,end-start]

def branchandbound(fulldata,powerbudget):
    start=perf_counter()
    N=len(fulldata)
    for i in range(len(fulldata)):
        if len(fulldata[i])==0:
            return [0,0]
    BestIntSoFar=greedy(fulldata,powerbudget,0)[0]
    Q=[]
    Q.insert(0,[0,0,0]) 
    while len(Q)>0:
        last=Q.pop(0)
        Values=fulldata[last[0]]
        for Pair in Values:
            Used=Pair[0]+last[1]
            Rate=Pair[1]+last[2]
            if Used>powerbudget:
                continue
            if last[0]==N-1:
                BestIntSoFar=max(BestIntSoFar,Rate)
            else:
                GreedyResult=greedy(fulldata,powerbudget-Used,last[0]+1)
                BestPossible=GreedyResult[1]
                if BestPossible+Rate>BestIntSoFar:
                    Q.insert(0,[last[0]+1,Used,Rate])
                BestIntSoFar=max(BestIntSoFar,Rate+GreedyResult[0])
    end=perf_counter()
    return [BestIntSoFar,end-start]                    
        
def online(N,M,K,pmax,rmax,powerbudget):
    Taken=[False for i in range(N)]
    fulldata=[[] for i in range(N) ] #We need this table to apply dynamic programming
    Rate=0
    usedpower=0
    for k in range(K):
        #We generate the random pairs
        Values=[[random.randint(1,pmax),random.randint(1,rmax),n] for e in range(M) for n in range (N)]
        #We add each pair to its respective channel in thetable fulldata
        for item in Values:
            Channel=item[2]
            fulldata[Channel].append([item[0],item[1]])
        #We sort the pairs by efficiency in a descending order
        Values=[[item[0],item[1],item[1]/item[0],item[2]]  for item in Values ]
        Values = sorted(Values, key=lambda x:-x[2])
        #We allocated a powerbudget and a maximum number of channels for each given user
        Channels=0
        #We iterate through the pairs until we have exceeded our powerbudget or maximum number of allocated channels
        for pair in Values:
            if not(Taken[pair[3]]) and  Channels<=N/K and usedpower+pair[0]<=powerbudget:
                Channels+=1
                Taken[pair[3]]=True
                Rate+=pair[1]
                usedpower+=pair[0]
    #We get the exact integer optimal solution using the dynamic programming algorithm
    DP=DP1(fulldata,powerbudget)
    return (Rate/DP[0],usedpower,DP[1])

def DP1(fulldata,powerbudget):
    preprocessq2(fulldata,powerbudget)
    preprocessq3(fulldata)
    preprocessq4(fulldata)
    return dynamic1(fulldata,powerbudget)
    
def main(path,testcode):
    print(path+":")
    
    #First we extract the table of powers, the table of rates and the powerbudget from the files
    data=parsefile(path)
    datap=data[0]
    datar=data[1]
    powerbudget=data[2]
    #We create one table of N arrays that contains pairs [Power,Rate] 
    fulldata=[[ [datap[i][e],datar[i][e]] for e in range(len(datap[i])) ] for i in range(len(datap))]
    unprocesseddata= deepcopy(fulldata)
    #We remove the impossible terms
    preprocessq2(fulldata,powerbudget)
    table=[len(x) for x in fulldata]
    print("Terms remaining after removing impossible terms: "+str(sum(table)))
    

    #We remove the IP-Dominated terms
    preprocessq3(fulldata)
    table=[len(x) for x in fulldata]
    print("Terms remaining after removing IP-dominated terms: "+str(sum(table)))
   
    #We remove the LP-dominated terms
    preprocessq4(fulldata)
    table=[len(x) for x in fulldata]
    print("Terms remaining after removing LP-dominated terms: "+str(sum(table)))
        
    if testcode==0 :
        
        
        #We apply the greedy algorithm on the preprocessed table
        greedyresult=greedy(fulldata,powerbudget,0)
        print("The rate achieved through the greedy algorithm is: "+str(greedyresult[1]))
        print("The algorithm used the power "+str(greedyresult[2])+" out of "+str(greedyresult[3]))
        print("The CPU runtime for the greedy algorithm is: " +str(greedyresult[4]))
        
        #We use LPsolver provided by scipy.optimize.linprog to compare with the greedy result
        lpsolved=LPsolveprocessed(fulldata,powerbudget)
        print("The optimal rate given by the LPsolver algorithm is: "+str(lpsolved[0]))
        print("The CPU runtime for the LPsolver algorithm is: "+str(lpsolved[1]))
    
    if testcode==1:    
        #We apply the first dynamic programming algorithm on the preprocessed table
        dynamicresult=dynamic1(fulldata,powerbudget)
        print("The optimal rate given by dynamic programming algorithm 1 is: "+str(dynamicresult[0]))
        print("The CPU runtime for the DP1 algorithm is: "+str(dynamicresult[2]))
    
    
        #We apply the second dynamic programming algorithm on the preprocessed table
        dynamicresult=dynamic2(fulldata,powerbudget)
        print("The optimal rate given by dynamic programming algorithm 2 is: "+str(dynamicresult[0]))
        print("The CPU runtime for the DP2 algorithm is: "+str(dynamicresult[1]))
        
        #We apply the Branch and Bound algorithm on the preprocessed table
        Bandb=branchandbound(fulldata,powerbudget)
        print("The optimal rate given by branch and bound algorithm is: "+str(Bandb[0]))
        print("The CPU runtime for the branch and bound algorithm is: "+str(Bandb[1]))
    
def onlinerun():
    #Let's run the offline algorithm 
    Off=[online(4,2,10,50,100,100) for i in range (10000)]
    Ratios=[item[0] for item in Off]
    UsedPowers=[item[1] for item in Off]
    UsedPowersDP=[item[2] for item in Off]
    print("The ratio achieved through the online algorithm is: "+str(mean(Ratios)))
    print("The online algorithm used an average power budget of: "+str(mean(UsedPowers)))
    print("The Corresponding offline algorithm used an average power budget of: "+str(mean(UsedPowersDP)))
    

    

paths=["test1.txt","test2.txt","test3.txt","test5.txt","test4.txt"]


