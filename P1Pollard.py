import numpy as np
import random
import math
from math import floor
from random import randint
import time 
#%%
# Test pour savoir si le nombre premier ou non 
def NombrePremier(p):
    if(p == 2):
        return True,1

    if(p == 1 or p%2==0):
        return False,2

    return Rabin_Miller(p, 50,2)

# Rabin Millier Test Nombre Premier
def Rabin_Miller(p, iterations,count):
    r=0
    s=p-1
    p=int(p)
    
    Premier=True
    count+=1
    while(s%2==0):
        count+=1 
        r+= 1
        s=int(s// 2)

    for i in range(0,iterations):
        a = randint(2,p-1)
        x = pow(a, s, p)
        
        count+=3
        if(x==1 or x==(p-1)):
            continue
        
        for j in range(r-1):
            x = pow(x,2,p)
            count+=1
            if(x==p-1):
                break
        else:
            Premier=False
            
    return Premier,count

# Algo Sqrt pour de très grands nombres, en utilisant la méthode de Newton
def large_sqrt(n): 
    x = n
    y = (x + 1) // 2
    count=0
    while(y < x):
        count+=1
        x = y
        y = (x + n // x) // 2
    count+=1
    return x,count   


def Erasthothene(n):
    count=0
    # créer une liste l de 2 à n 
    L=list(range(2,n+1))
    count+=1
    
    i=2
    
    while(i<=large_sqrt(n)[0]):
       count+=1
       # Si i dans la liste, sinn supprimer ses multiples
       count+=1
       if(i in L):
           #j donne des multiples de i, on commence de 2*i et on incrémente i par i (exemple i=2 on fait +2)
           for j in range(i*2,n+1,i):
               #Suppression des multiples s'ils sont trouvés dans la liste
               count+=1
               if(j in L):
                   count+=1
                   L.remove(j)
       i=i+1
    count+=1
    return L,count
#%%

#----------------------------------------------------------------------
#-------                    ALGO P-1 POLLARD                   --------    
#----------------------------------------------------------------------   

list_fact=[]

def p_1_Pollard(n,B,count):
   
    # pas de diviseur premier 1
    if (n == 1): 
        count+=1 # if (n == 1): 
        #print("No Prime divisor for 1")
        list_fact.append(n)
        count+=1 # list_fact.append(n)
        return list_fact,count
    
    if(NombrePremier(n)[0]==True):
        list_fact.append(n)
        count+=1
        return list_fact,count
    
    # 1) Choisir Seuil de Friabilité (dans les paramètres)
    
    BoundOk=False
    RandomOk=False
    
    while(BoundOk==False):
        count+=1 #  while(BoundOk==False):
        
        # 1)
        #print("Bound value B : "+str(B))
        
        # Reinitialisation
        RandomOk=False
            
        ## Condition d'arrêt B < sqrt(N)
        count+=1 #  if(B> (np.sqrt(n)*100)):
        if(B> (large_sqrt(n)[0]*100)):
            #print("Algorithme terminé ! B > sqrt(n)")
            return -1
       
        if(B<2):
            count+=1
            #print("B Trop bas, Erreur !")
            return list_fact,count
        

        # 2) define M = product(q^floor( log(B)/log(q)) )
        M=1
        
        
        listE=Erasthothene(B)[0]
        
        for q in listE:
            M= M*(q**(floor(math.log(B)//math.log(q))) )
            #print("itérations "+str(q))
            #print("M : "+str(M))
        #print(M)
        
        """
        for q in range(2,B+1):
            if(NombrePremier(q)[1]==True):
                M= M*(q**(floor(math.log(B)//math.log(q))) )
                #print("itérations "+str(q))
                #print("M : "+str(M))
        #print(M)
        """
        
        while(RandomOk==False):           
            count+=1 #  while(BoundOk==False):
            
            # 3) randomly pick a coprime to n (fix a=2 if n is odd (impair))
         
            a=random.randint(2,n-1) # random number between 1 & n
            #print("Random number a is : "+str(a))
                
            # 4) compute g = gcd(a^M − 1, n) 
            g=np.gcd( (a**M) - 1,n)
            
            #print("gcd is : "+str(g))
            #print("")
                    
            # 5) 
            count+=1+1 #  if(g!=1 and g!=n):
            if(g!=1 and g!=n):
                RandomOk=True
                BoundOk=True
                #print("Facteur Non Trivial : "+str(g))
                count+=NombrePremier(g)[1]+1 # if(NombrePremier(g)[0]==True):
                if(NombrePremier(g)[0]==True):
                    list_fact.append(g)
                    count+=1 # list_fact.append(g)
                    #print("Current List Fact is : "+str(list_fact))
                    
                else:
                    other_factor=n//g
                    count+=NombrePremier(other_factor)[1]+1 # if(NombrePremier(other_factor)[0]==True):
                    if(NombrePremier(other_factor)[0]==True):
                        list_fact.append(other_factor)
                        count+=1 # list_fact.append(other_factor)
                        #print("Current List Fact is : "+str(list_fact))
                    else:
                        return p_1_Pollard(other_factor,B,count)
                   
                count+=1+1+1+NombrePremier(n//g)[1] # if(n%g==0  and NombrePremier(n//g)[0]==True):
                if(n%g==0 and NombrePremier(n//g)[0]==True):
                    list_fact.append(n//g)
                    count+=1 # list_fact.append(n//g)
                    #print("Current List Fact is : "+str(list_fact))
                
                # Facteur restant par exemple 11*11*13 on a fait 11*143 il faut verifier que 143 est encore divisible
                elif(n%g==0 and NombrePremier(n//g)[0]!=True):
                    count+=1+1+1+NombrePremier(n//g)[1] # if(n%g==0  and NombrePremier(n//g)[0]==True):
                    #print("Reccurence !")
                    #print("")
                    return p_1_Pollard(n//g,B,count)
                
                #print("Les facteurs sont : "+str(list_fact))

                return list_fact,count
                        
            count+=1 #  if(g==1):
            if(g==1):
                
                # Choisir un B plus grand
                count+=1 #  if(B<=np.sqrt(n)*pow(10,len(str(n)))):
                if(B<=large_sqrt(n)[0]*pow(10,len(str(n)))):
                    RandomOk=True
                    B=B+(len(str(n))//pow(10,len(str(n))))  # len(str(n)) longueur du chiffre n 
                                                        # Par exemple 12405 : len(n) = 5
                
                # Error
                else:
                    #print("Erreur g==1 !")
                    return list_fact,count
            
            count+=1 # if(g==n):
            if(g==n):
                
                count+=NombrePremier(g)[1]+1 # if(NombrePremier(g)[0]==True):
                if(NombrePremier(g)[0]==True):
                    list_fact.append(g)
                    count+=1 # list_fact.append(g)
                    #print("Listes des facteurs actuels: "+str(list_fact))
                
                # Error
                else:
                    if(B<=large_sqrt(n)[0]*pow(10,len(str(n)))):
                         RandomOk=True
                         B=B+(len(str(n))//pow(10,len(str(n))))
                    else:  
                        print("Failure ! ")
                        
                    #print("Les Facteurs sont : "+str(n))
                #return list_fact,count

#%%
def P1Pollard_Chrono(n,B):
    time1=time.time() 
    resultat=p_1_Pollard(n,B,0)
    time2=time.time() 
    time3=time2-time1
    return resultat,time3

n=146182562237
n=264839967043414254127
#n=2*2*2*5
#n=4217*4217*4217
#n=290919706205487829572593

B=20


B=5
#n=703253 
n=101*137

print("")
print("---------------")
print("- P-1 Pollard -")
print("---------------")
print("")
print("Nombre à Factoriser n = "+str(n))
print("")
result=P1Pollard_Chrono(n,B)
print("Liste des facteurs : "+str(result[0][0]))
print("Nombres d'Opérations : "+str(result[0][1]))
print("Temps : "+str(result[1]))





#%%
"""

def plot_P1(start,end,B):
    
    yop_p1=[]
    ychrono_p1=[]
    
    for i in range(start,end,10000):
        print("Iteration : "+str(i))
        result=P1Pollard_Chrono(i,B)
        print(result)
        
        yop_p1.append(result[0][1])
        ychrono_p1.append(result[1])
        print("")
        
    return [yop_p1,ychrono_p1]


import matplotlib as plt
import matplotlib.pyplot as plt


start=5
end=1000005

res_p=plot_P1(start,end,B)
yop_p=res_p[0]
ychrono_p=res_p[1]
      
x= [i for i in range(start,end,10000)]


plt.plot(x, yop_p,label = "P-1 Pollard")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Opérations nécessaires')
axes.set_title('P-1 Pollard')
plt.show()

plt.plot(x, ychrono_p,label = "P-1 Pollard")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Temps en secondes')
axes.set_title('P-1 Pollard')
plt.show()

"""
