import numpy as np
from math import floor
import numpy as np
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
#%%

liste_ANaif=[]

def Algo_Naif_Rec(nb,count):
    res=NombrePremier(nb)
    count+=res[1] # nombrePremier(n)
    if(res[0]):        
        liste_ANaif.append(nb)
        count+=1  # if(res[0]):   
        return liste_ANaif,count
    else:
        count+=1  # if(res[0]):   
        for k in range(2,floor(large_sqrt(nb)[0])+1):
        		
            res=NombrePremier(k)
            count+=1+res[1]  # if(res[0]):  + nombrePremier(n)
            if(res[0]):          
                # Si il est divisible, on ajoute k nombre premier dans 
                # la liste de la decomposition
                count+=1 # n%k==0
                if(nb%k==0):
                    liste_ANaif.append(k)
                    count+=1 # liste_ANaif.append(k)
                    return Algo_Naif_Rec(nb//k,count)
        return liste_ANaif,count

#%%
def AlgoNaif_Chrono(n):
    time1=time.time() 
    resultat=Algo_Naif_Rec(n,0)
    time2=time.time() 
    time3=time2-time1
    return resultat,time3

nbPremier=2*3*5*7*11*2*2*5

n=146182562237
n=1000
n=264839967043414254127

n=10153331*373587883*76695841

n=703253 

n=3480
print("")
print("-------------")
print("- Algo Naïf -")
print("-------------")
print("")
print("Nombre à Factoriser n = "+str(n))
print("")
result=AlgoNaif_Chrono(n)
print("Liste des facteurs : "+str(result[0][0]))
print("Nombres d'Opérations : "+str(result[0][1]))
print("Temps : "+str(result[1]))

#%%
"""
def plot_naif(start,end):
    
    yop_rho=[]
    ychrono_rho=[]
    
    for i in range(start,end,10000):
        print("Iteration : "+str(i))
        result=AlgoNaif_Chrono(i)
        print(result)
        
        yop_rho.append(result[0][1])
        ychrono_rho.append(result[1])
        print("")
        
    return [yop_rho,ychrono_rho]


import matplotlib as plt
import matplotlib.pyplot as plt


start=5
end=1000005

res_naif=plot_naif(start,end)
yop_naif=res_naif[0]
ychrono_naif=res_naif[1]
      
x= [i for i in range(start,end,10000)]


plt.plot(x, yop_naif,label = "Naif")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Opérations nécessaires')
axes.set_title('Algorithme Naif')
plt.show()

plt.plot(x, ychrono_naif,label = "Naif")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Temps en secondes')
axes.set_title('Algorithme Naif')
plt.show()


"""




"""
print("")
print("--------------------------------------------------------")
print("--                 DEBUT ALGORITHM NAIF                 --")
print("--------------------------------------------------------")
print("")
print("")
print("")
print("Test Algo Naif")
print("Valeur t : "+str(nbPremier))
resultat=Algo_Naif_Rec(nbPremier,0)
print("Decomposition en Facteur Premier : "+str(resultat[0]))
print("Compteur : "+str(resultat[1]))
print("")
print("--------------------------------------------------------")
print("--                 FIN ALGORITHM NAIF                 --")
print("--------------------------------------------------------")
print("")
print("")
"""
