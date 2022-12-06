import random as rnd
import numpy as np


# Randomize function with a restriction in the order
def resrand(vec,maxseq):
    
    #if maxseq = 3, it check for 4 numbers in an order
    maxseq += 1

    #shuffle untill sequence does not contain more tha maxseq in a row
    passed = False
    while passed == False:
        
        #new assignment each time to get track of order
        newvec=vec[:]
        
        #define as true
        passed = True
    
        #randomize
        r_ind=rnd.sample(range(0, len(newvec)), len(newvec))
        newvec = [newvec[x-1] for x in r_ind]
        
        
        #loop over sequences
        for i in range(maxseq,len(newvec)+1):

            #check if a sequence diff adds up to 0 for example > [4,4,4,4]                    
            if sum(abs(np.diff(newvec[i-maxseq:i])))==0:
                
                #define false and loop continues
                passed = False
    return newvec,r_ind

def makeiti(stimcat,itiminmax):
    
    l_iti = [0] * len(stimcat)
    
    #settings
    v_min = itiminmax[0]
    v_max = itiminmax[1]
    
    #for each category make a list
    v_unique = list(set(stimcat))
    
    #loop over categories    
    for c_cat in v_unique:
        
        #create a fixed list
        temp_l_iti = np.linspace(v_min,v_max,stimcat.count(c_cat))
        temp_l_iti=[int(x) for x in temp_l_iti]
        
        #shuffle for each category
        rnd.shuffle(temp_l_iti)
        
        #get position in orginal list
        pos = [x for x, y in enumerate(stimcat) if y == c_cat]
        
        #loop over itis
        for c_pos in range(len(pos)):
            
            #place iti in final list
            l_iti[pos[c_pos]] = temp_l_iti[c_pos]

    #return randomized ITI list with equal ITIs for all categories
    return l_iti

#print a list to a txt file
def savetxt(thelist,filename):
    
    #thelist=list(map(list, zip(*thelist)))
    
    thefile = open(filename, 'w')
    
    for item in thelist:
        thefile.write("%s\n" % item)
    
    thefile.close()

#read in txt file 
def readtriallist(filename):


    triallist=[ [], [], [], [], [], [], [], [], [], []]
    with open(filename) as f:
        for line in f:
            dat=line.replace("'","").replace("[","").replace("]","").replace("\n","").split(', ')
            for i,val in enumerate(dat):
                triallist[i].append(val)       
    
    return triallist 



def divide_chunks(l, n): 
      
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 