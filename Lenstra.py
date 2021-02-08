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
#%%
"""
    Factorisation de Lenstra par les courbes elliptiques
    
    Algo Probabiliste rapide pour la décomposition en produit de facteurs premiers
    qui emploie les courbes elliptiques
    
    Amélioration de l'algo p-1 Pollard
    
    
    Etapes de L'Algorithme
    
    1)
    Prendre une courbe elliptique aléatoire sur Z avec un point A sur elle. 
    Alors, nous considérons la loi de groupe sur cette courbe modulo n - ceci est possible 
    car la plupart des résidus modulo n ont des inverses, 
    qui peuvent être trouvés en utilisant l'algorithme d'Euclide 
    et en trouvant un résidu non-inversible équivalent à la factorisation de n
    
    2) 
    Calculer eA dans ce groupe, où e est le produit de petits nombres premiers 
    élevés aux petites puissances, comme dans la méthode p−1 de Pollard. 
    Il peut donner un nombre premier en une fois, et est ainsi efficace.
    
    3)
    Avec un peu de chance, eA est l'élément nul du groupe de la courbe elliptique dans Fp, 
    mais pas dans Fq pour un autre diviseur premier q de n 
    (comme dans la méthode p−1 de Pollard, il est très improbable que les deux groupes 
    aient un ordre qui soit un diviseur de e). 
    Alors nous pouvons trouver un facteur de n 
    en calculant le PGCD de la première coordonnée de A et n, 
    car cette coordonnée sera nulle dans Fp
    
    4)
    Si cela ne marche pas, 
    il suffit de recommencer avec une autre courbe ou un autre point de départ.


    Complexité
    
    Dépend de la taille du facteur.
    Elle peut être exprimée par O(e(√2 + o(1)) √ln p ln ln p)
    où p est le plus petit facteur de n.

"""
#%%

# Determine l'inverse d'un nombre dans Z/nZ / Algo + optimisé
# Algorithme Euclide Etendu  
def extended_gcd(a, b):
    count=0
    r, reste = abs(a), abs(b)
    x, u, y, v = 0, 1, 1, 0
    while reste:
        count+=1 # while reste:
        r, (quotient, reste) = reste, divmod(r, reste)
        x, u = u - quotient*x, x
        y, v = v - quotient*y, y
    return r, u * (-1 if a < 0 else 1), v * (-1 if b < 0 else 1),count


# P=(x1,y1)
# Q=(x2,y2)

# Addition de points d'une courbe elliptique dans Z/nZ
def Addition_Point(P,Q,Fp):
    count=0
    count+=1 # if(P =='O'): 
    if(P =='O'): 
         return Q,count
    count+=1 # if(Q =='O'):
    if(Q =='O'): 
        return P,count
    
    # x1 = x2 && y1 = -y2 => P = -Q
    count+=1+1+1+1+1+1+1 # if ( P[0]%Fp == Q[0]%Fp) and ( P[1]%Fp == -Q[1]%Fp): 
    if ( P[0]%Fp == Q[0]%Fp) and ( P[1]%Fp == -Q[1]%Fp): 
        return 'O',count
    
    # Si x1 = x2 && y1 = y2  => P = Q 
    count+=1+1+1+1+1+1+1 # if (P[0]%Fp == Q[0]%Fp) and ( P[1]%Fp == Q[1]%Fp):    
    if (P[0]%Fp == Q[0]%Fp) and ( P[1]%Fp == Q[1]%Fp):       
        
        a=randint(2,Fp-1)
        
        # m= (3*x1² + a) / 2y1
        num= ( (3*(P[0]**2) )+ a )%Fp
       
        ### Renvoie Inverse de denom ou retourne denom si il n'a pas d'inverse (utile pr Lenstra)
        
        # Euclide Etendu 
        denom=2*P[1]
        count+=1 #  P[1]
        g,x,y,c=extended_gcd(denom,Fp)
        count+=c # extended_gcd(denom,Fp)
        count+=1 #if(g!=1):
        if(g!=1):
            #print("Pas d inverse pour : "+str(denom)+"%"+str(Fp))
            return [denom,'No Inverse'],count ### Lenstra retourne d 
        
        #On inverse denom si il a un inverse
        denom_inv=(x%Fp)
        
        m=num*denom_inv
        
        # x3 = m²-2x1
        x3=( (m**2)-2*P[0]) %Fp  
        count+=1 # P[0]
              
    # Si x1 != x2 && y1!=y2 => P != Q  
    else :
        # m=(y2-y1)/(x2-x1)
        num=( Q[1] - P[1])%Fp
        count+=1+1 #  Q[1]  P[1]
        
        ### Renvoie Inverse de denom ou retourne denom si il n'a pas d'inverse (utile pr Lenstra)      
        # Euclide Etendu 
        count+=1+1 # Q[0] - P[0]
        denom=Q[0] - P[0]
        g,x,y,c=extended_gcd(denom,Fp)
        count+=c # extended_gcd(denom,Fp)
        count+=1 #if(g!=1):
        if(g!=1):
            #print("Pas d inverse pour : "+str(denom)+"%"+str(Fp))
            return [denom,'No Inverse'],count
  
        #On inverse denom si il a un inverse      
        denom_inv=(x%Fp)
    
        m=num*denom_inv
        
        # x3 = m²-x1-x2
        count+=1+1 # P[0]-Q[0] ) 
        x3 =( (m**2)-P[0]-Q[0] ) 

    # (x3, -y1 + m*(x1-x3) )
    count+=1+1 #  P[0] - x3 ) - P[1])
    return ([ x3%Fp ,( m *( P[0] - x3 ) - P[1]) %Fp]),count


# Test For Addition 
"""
print("Test Addition ")

a=2
Fp=97

P=[15,6]
Q=[5,9]

print("Coord use : P = "+str(P)+" , Q = "+str(Q)+" , n = "+str(Fp) )
print(Addition_Point(P,Q,Fp) )
print("")


# No inversible
P=[2,11]
Q=[2,11]
Fp=36
print("No Inversible Test for lenstra")
print("Coord use : P = "+str(P)+" , Q = "+str(Q)+" , n = "+str(Fp) )
print(Addition_Point(P,Q,Fp) )
print("")

print("Test Inversible")
print(modinv(3,Fp))
print("")
"""

#----------------------------------------------------------------

#----------------------------------------------------------------


#%%
# Multiplication nP (Addition nP = P + P + P + ... + P n fois)
def Multiplication_Points(n,P,Fp):
    result = 'O'
    pow_2P = P
    count=0
    while(n!=0):     
        count+=1 # while(n !=0):
        count+=1+1 # if(pow_2P[1]!='No Inverse'):
        
        if(pow_2P[1]!='No Inverse'):
            if (n%2==1):
                    res = Addition_Point( pow_2P , result,Fp)

                    result=res[0]

                    count+=res[1]
                    count+=1+1 #if(result[1]=='No Inverse'):
                    if(result[1]=='No Inverse'):
                        return result,count
            res=Addition_Point( pow_2P , pow_2P,Fp)
            pow_2P =res[0]
            count+=res[1]
            #print("n = "+str(n))
            #print("pow_2P = "+str(pow_2P))
            #print("result = "+str(result))
            #print("")
            n=n//2
        else:
            return pow_2P,count
            
    return result,count

"""
print("Test Mult No inverse")
print(Multiplication_Points(2,P,Fp))
print("")
"""

"""
print("Test Multiplication nP")
print("")

a=2
b=3
n=4
Fp=1999
P=[3,6]

print("Values are : n = "+str(n)+" , a = "+str(a)+" , Fp = "+str(Fp)+" , P = "+str(P))
print("nP value is : "+str(Multiplication_Points(n,P,Fp) ) )

"""


# Algorithme Crible Erasthothène  
# Retourne la liste de tous les nombres premiers <= n
def Erasthothene(n):
    count=0
    # créer une liste l de 2 à n 
    L=list(range(2,n+1))
    count+=1
    
    i=2
    
    while(i<=np.sqrt(n)):
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

lenstra_list_factor=[]
#----------------------------------------------------------------------
#-------                    ALGO LENSTRA                   ------------    
#---------------------------------------------------------------------- 
def Lenstra(n,B,count):
    
     # 1) Vérifier que :
     # n est pas premier
     count+=1
     if(NombrePremier(n)[0]==True):
         return n,count
     
     # n n'est pas divisible par 2
     count+=1
     if(n%2==0):
         return 2,count
     
     # n n'est pas divisible par 3
     count+=1
     if(n%3==0):
         return 3,count
     
        
     if(B<2):
         #print("B : "+str(B))
         print("Algo Fail ! ")
         return lenstra_list_factor,count
        
        
     # 2) Choisir des entiers aléatoires a,x,y entre 1 et n
     RandomOk=False
     while(RandomOk==False):
         count+=1
         Steptwo=False
         x,y=randint(1,n-1),randint(1,n-1)
         
         while(Steptwo==False):
             count+=1
             
             a=randint(1,n-1)
             print("a : "+str(a))
             print("x : "+str(x))
             print("y : "+str(y))
             
             
             
             # 3) Calculer b = y**2 − x**3 − a*x (mod n).
             b= (pow(y,2)-pow(x,3)-a*x) %n
             print("b : "+str(b))
             print("")
             # 4) Calculer d= pgcd(4a**3 + 27b**2, n)
             f=int(4*pow(a,3) + 27*pow(b,2))%n
             #print(f)
             d=np.gcd(f,n)
             print("d : "+str(d))
             count+=1+1
             # Si 1<d<n, retourner le facteur non trivial
             if(1<d and d<n):
                print("d : "+str(d))
                count+=1+NombrePremier(int(d))[1] # if(NombrePremier(int(facteur))==True):
                if(NombrePremier(int(d))[0]==True):
                    lenstra_list_factor.append(int(d))
                    count+=1 # lenstra_list_factor.append(int(facteur))
                
                    # tous les facteurs sont premiers on a finis
                    count+=1+NombrePremier(int(n//d))[1] # if(NombrePremier(int(facteur))==True):
                    if(NombrePremier(int(n//d))[0]==True):
                        lenstra_list_factor.append(int(n//d))
                        count+=1 # lenstra_list_factor.append(int(n//facteur))
                        print("Listes des facteurs non trivials :")
                        return lenstra_list_factor,count
                
                # il reste encore à décomposer en nombre premier
                    else:
                        return Lenstra(int(n//d),B-5,count)
                               
                    # Le facteur trouvé n'est pas premier
                else:
                    autre_facteur=int(n//d)
                    count+=1+NombrePremier(autre_facteur)[1] # if(NombrePremier(autre_facteur)[0]==True):
                    if(NombrePremier(autre_facteur)[0]==True):
                    
                       lenstra_list_factor.append(autre_facteur)
                       count+=1 # lenstra_list_factor.append(autre_facteur)
                       return Lenstra(int(n//autre_facteur),B-5,count)
                
                    else:
                        return (Lenstra(int(np.gcd(d,n)),B-5,count))
            
            
             # Si d=1, aller à l'étape 5 
             count+=1
             if(d==1):
                 Steptwo=True
                 
             # Si d=n, aller à l'étape 2, choisir un a différent
             
             count+=1
             if(d==n):
                 #print("Retour à l étape 2")
                 print("Choisir un a différent")
                
                 
         # 5) Courbe Elliptic : y**2 = x**3 + ax + b & P=[x,y]
         P=[x,y]
         print("P : "+str(P))
         # 6) Choisir k = lcm(2,3,...B)
         k=1
         liste_premier=Erasthothene(B)[0]
         #print("Liste Nombre premiers : "+str(liste_premier))
         for p in range(0,len(liste_premier)):
             count+=1
             k=np.lcm(k,liste_premier[p])
             
         print("k : "+str(k))
         print("")
         # 7) Compute kP % n
         res=Multiplication_Points(k,P,n)
         P=res[0]
         count+=res[1]
         print("kP : "+str(P))
         
         count+=1
         if(P[1]=='No Inverse'):
            d=P[0]
            print("No inverse d : "+str(d))
            
            count+=1
            if(np.gcd(d,n)<n): 
                
                facteur=np.gcd(d,n)
                count+=1+NombrePremier(int(facteur))[1] # if(NombrePremier(int(facteur))==True):
                if(NombrePremier(int(facteur))[0]==True):
                    lenstra_list_factor.append(int(facteur))
                    count+=1 # lenstra_list_factor.append(int(facteur))
                
                    # tous les facteurs sont premiers on a finis
                    count+=1+NombrePremier(int(n//facteur))[1] # if(NombrePremier(int(facteur))==True):
                    if(NombrePremier(int(n//facteur))[0]==True):
                        lenstra_list_factor.append(int(n//facteur))
                        count+=1 # lenstra_list_factor.append(int(n//facteur))
                        
                        #print("Listes des facteurs non trivials :")
                        print("")
                        print("-------------------------------------")
                        print("")
                        return lenstra_list_factor,count
                
                # il reste encore à décomposer en nombre premier
                    else:
                        return Lenstra(int(n//facteur),B-13,count)
                               
                    # Le facteur trouvé n'est pas premier
                else:
                    autre_facteur=int(n//facteur)
                    count+=1+NombrePremier(autre_facteur)[1] # if(NombrePremier(autre_facteur)[0]==True):
                    if(NombrePremier(autre_facteur)[0]==True):
                    
                       lenstra_list_factor.append(autre_facteur)
                       count+=1 # lenstra_list_factor.append(autre_facteur)
                       return Lenstra(int(n//autre_facteur),B-5,count)
                
                    else:
                        return (Lenstra(int(np.gcd(d,n)),B-5,count))
                    
            # Found d factor of n
            #return np.gcd(d,n)
            
            count+=1 # if(np.gcd(d,n)==n):        
            if(np.gcd(d,n)==n):
                #print("Réduire k !")
                k=k-len(str(k))
                #k=k-1
         else:
             #print("Choisir une valeur différente pour a,x,y")
             Steptwo=True 
#%%
def Lenstra_Chrono(n,B):
    time1=time.time() 
    resultat=Lenstra(n,B,0)
    time2=time.time() 
    time3=time2-time1
    return resultat,time3

n=146182562237
#n=1811706971
#n=5*7*11
B=47


#n=290919706205487829572593
n=264839967043414254127

n=703253 
B=25

n=455839
B=5
n=1081

print("")
print("Lenstra")
print("")
print("Factorisation de "+str(n))   
print(str(Lenstra_Chrono(n,B)))




print("")
print("-----------")
print("- Lenstra -")
print("-----------")
print("")
print("Nombre à Factoriser n = "+str(n))
print("")
result=Lenstra_Chrono(n,B)
print("Liste des facteurs : "+str(result[0][0]))
print("Nombres d'Opérations : "+str(result[0][1]))
print("Temps : "+str(result[1]))


"""
print("Test Algo Lenstra")
print("n : "+str(n))
print("Limite B : "+str(B))
resultat=Lenstra_Fct(n,B,0)
print("Decomposition en Facteur Premier : "+str(resultat[0]))
print("Compteur : "+str(resultat[1]))
"""


#%%
"""
B_len=30
def plot_len(start,end,B):
    
    yop_p1=[]
    ychrono_p1=[]
    
    for i in range(start,end,10000):
        print("Iteration : "+str(i))
        result=Lenstra_Chrono(i,B)
        print(result)
        
        yop_p1.append(result[0][1])
        ychrono_p1.append(result[1])
        print("")
        
    return [yop_p1,ychrono_p1]


import matplotlib as plt
import matplotlib.pyplot as plt


start=5
end=1000005

res_len=plot_len(start,end,B)
yop_len=res_len[0]
ychrono_len=res_len[1]
      
x= [i for i in range(start,end,10000)]


plt.plot(x, yop_len,label = "Lenstra")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Opérations nécessaires')
axes.set_title('Lenstra')
plt.show()

plt.plot(x, ychrono_len,label = "Lenstra")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Temps en secondes')
axes.set_title('Lenstra')
plt.show()

"""


"""
print("")
print("--------------------------------------------------------")
print("--               ALGORITHME LENSTRA                   --")
print("--------------------------------------------------------")
print("")


#n=13*17*2*19
n=146182562237
#n=34665061


print("")
print("--------------------------------------------------------")
print("--             FIN ALGORITHME LENSTRA                 --")
print("--------------------------------------------------------")
print("")
"""