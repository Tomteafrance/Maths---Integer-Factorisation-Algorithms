import numpy as np
import math
from itertools import chain
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
##################################################################
## -------------------  INITIALISATION   -------------------------
##################################################################    

def Seuil_Friabilite(n):
    B=int((np.exp((1/2)*large_sqrt(math.log(n)*math.log(math.log(n)))[0])))
    return B

def Exponentiation_rapide(x,n):
    count_Legendre.append(1)
    if(n==1):
        
        return x
    if(n%2==0):
        return Exponentiation_rapide(x**2,n//2)
    count_Legendre.append(1)
    if(n%2==1):
        count_Legendre.append(1)
        return x*Exponentiation_rapide(x**2,(n-1)//2)


def Legendre(a,p):
    return Exponentiation_rapide( a,(p-1)//2 )%p

# Réduction de la Base de facteurs, on utilise seulement les nombres premiers qui sont des 
# résidue quadratiques
def Base_Facteurs(N,B):
    count=0
    
     # let t = TT(B)  liste de nombres premiers en dessous de B 
    liste_premier=Erasthothene(B)
    t=liste_premier[0]
    count+=liste_premier[1]
    base_facteur=[]
    base_facteur.append(2)
    count+=1
    
     # Restriction base des facteurs 
    for i in range(1,len(t),1):
         if(Legendre(N,t[i])==1):
             base_facteur.append(t[i])
             count+=1
    return base_facteur,count
#%%
##################################################################
## -------------------      CRIBLAGE     -------------------------
################################################################## 

# Algorithme de Shanks Tonelli : Trouve les racines des équations du type r²=N %p
def Shanks_Tonelli(n, p):
    count=0
    assert Legendre(n, p) == 1, "Pas un résidu quadratique !"  
    # Si la condition n'est pas vrai, print(Pas un résidu quadratique)
    Q = p-1
    S = 0
    
    while(Q%2 == 0):
        count+=1
        Q = Q// 2
        S= S+1
    count+=2
    if(S==1):
        R = pow(int(n), (p + 1) // 4, p)
        return [R,p-R],count
    for z in range(2, p):
        if(p-1 == Legendre(z, p)):
            break
    M = S
    #c = pow(z, q, p)   # c=z**q % p
    c=(z**Q)%p
    #t = pow(n, q, p)
    t=(n**Q)%p
    #r = pow(n, (q + 1) // 2, p)
    R= (n**((Q+1)//2)) %p
    
    d = 0
    while( (t-1) % p!=0):
        count+=1
        d = (t**2) % p
        for i in range(1, M):
            if( (d-1)%p == 0):
                break
            d = (d**2) % p
            
        if(M-i-1<0):
            cpower=1
        else:
            cpower=M-i-1
            
        b= pow(c,pow(2,cpower,p),p)
        #b = pow(c, int(np.abs(1 << (M-i-1))), p) # b=c**2M-i-1
        M = i

        c = (b **2) % p          # c= b² %p
        t = (t * b**2) % p       # t= tb² % p
        R = (R * b) % p          # R=Rb%p
    return [R,p-R],count

# Crible pour les congruence, génère une liste d'équation x²-N en partant de sqrt(N)
def liste_crible(N,Intervalle):
    count=0
    liste_congruence=[]
    
    racine_N=large_sqrt(N)[0]
    
    for i in range(racine_N,racine_N+Intervalle):
        liste_congruence.append(i**2-N)
        count+=1
    return liste_congruence,count

# Recherche des congruences a² congru b² % n, I = B**3
def Criblage(base_facteur,N,Intervalle):
    count=0
    
    racine_N=large_sqrt(N)[0]
    
    # On redéfinit la séquence de criblage        
    liste_congruence=liste_crible(N,Intervalle)[0]
    count+=liste_crible(N,Intervalle)[1]
        
    # On garde une copie la liste cribler pour plus tard
    liste_cribler = liste_congruence.copy() 
    print("liste congruence : "+str(liste_cribler))
      
    # Si le premier facteur premier est 2 dans la base de facteur  
    # Cas particulier pour le facteur 2 
    count+=2
    if(base_facteur[0] == 2):
        i = 0
        # Tant que ce n'est pas divisible  on continue à tester tous les cribles jusqu'à 
        # qu'on trouve le premier terme 
        while(liste_cribler[i]%2!=0):
            count+=2
            i=i+1
        count+=2 
        # de i à la longueur de la liste criblage,
        # les pas sont de 2 car on a par exemple x=46 on veut aller à 48 50 puis 52 ...
        
        for j in range(i,len(liste_cribler),2): 
            # Tant que c'est divisible par 2, on divise encore et encore
            
            while(liste_cribler[j]%2==0):
                count+=2
                liste_cribler[j] =liste_cribler[j]// 2
                count+=2
            count+=2
            
    # Cas Général, Congruence => Algorithme de Shanks Tonelli
    # On parcourt la liste de la base de facteurs premiers (on saute 2)
    
    for premier in base_facteur[1:]: 
        count+=1
        solutions_residues = Shanks_Tonelli(N,premier)[0]
        count+=Shanks_Tonelli(N,premier)[1]
        # cela nous donne les solutions à l'équations x²=n % p où p est le facteur premier 
      
        print("Solutions résidues quad "+str(premier)+" : "+str(Shanks_Tonelli(N,premier)[0]))
        # On parcourt les solutions 
        for r in solutions_residues:
            count+=1
            print("premier : "+str(premier))
            print("r : "+str(r))
            # Parcourt solution - racine de n % p à la longueur liste crible, 
            # pas : facteur premier de la base
            for k in range((r-racine_N)%premier, len(liste_cribler),premier):
                count+=1
                #print("k : "+str(k))
                # Tant que que c'est divisible pour tous les pas p nombres 
                while(liste_cribler[k]%premier==0):
                    count+=2
                    liste_cribler[k]=liste_cribler[k]// premier    
                
                count+=2

    liste_num_friables = []  # liste des nombres x²-N 
    a_list = [] #liste a_1, a_2,...a_n les nombres d'origines en partant de la racine par exemple racine n  = 45, ce sera 46 47 ....
    index = []  # indices de ces nombres 
    
    print("")
    print("liste cribler : "+str(liste_cribler))
    print("")
    for i in range(len(liste_cribler)): 
        count+=1
        # Pas de solutions 
        #probabilité Pas de solutions is 2^-T, T facteur de tolérance IMPORTANT
        count+=1
        if len(liste_num_friables) >= len(base_facteur)+T: 
            break
        
        # Si on trouve 1 nombre lisse  [on a que des 1 ou 0 normalement apres toutes les divisions]
        count+=1
        if(liste_cribler[i] == 1): # On trouve un nombre B-friable
            liste_num_friables.append(liste_congruence[i]) # on ajoute x dans la liste des nombres friables
            count+=1
                                             
            a_list.append(i+racine_N) # on ajoute i + racine_N représente le nombre originel, 
            count+=1            
            index.append(i)                # on ajoute les indices de ce nombre
            count+=1

    return(liste_num_friables,a_list,index,count)

#%%
##################################################################
## ------------------  ALGEBRE LINEAIRE   ------------------------
##################################################################    

#Liste des facteurs premiers d'un nombre n, par exemple 46=2*23 => return [2,23]
def decompo_facteur(x,base_facteur):
    liste_facteur=[]
    count=1
    if(x<0):
        liste_facteur.append(-1)
        count+=1
    for premier in base_facteur:
        count+=1
        if(premier==-1):
            pass
        else:
            while(x%premier==0):
                count+=1
                liste_facteur.append(premier) # On rajoute le facteur par exemple 46%2 46= 2*..*. puis on a 23
                count+=1
                x//=premier
            count+=1
    return liste_facteur,count


def Matrice(liste_num_friables,base_facteur):
# Génère les vecteurs des exposants % 2, de la liste des nombres lisses qu'on a obtenu
# Puis créer la matrice
    
# Vecteur dont les composantes sont les puissances des nombres premiers 
# L'idée c'est que pour chaque facteur premier dans la base de facteur
# On regarde combien il y a de 2,3,5,7... len(base_facteur) dans a_1,a_2....,a_n
# puis on la matrice M se composera de 
# ligne 1 : les exposants de 2 pour chaque nombre 
# ligne 2: les exposants de 3 etc...
# Puis on transposera la matrice pour avoir chaque vecteur en colonne
    count=0
    M = []
    carre=False
    base_facteur.insert(0,-1)
    for x in liste_num_friables:
        print("Nombre friable : "+str(x))
        
        
        liste_exposant = [0]*(len(base_facteur)) # initialise les exposants des exposants du vecteur 
        x_liste_decompo = decompo_facteur(x,base_facteur)[0] # les facteurs de n 
        count = decompo_facteur(x,base_facteur)[1]
        print("decomposition facteur : "+str(x_liste_decompo))
        for i in range(len(base_facteur)):
            count+=2
            if(base_facteur[i] in x_liste_decompo):
                count+=3
                liste_exposant[i] = (liste_exposant[i] + x_liste_decompo.count(base_facteur[i])) % 2
                # On fait regle les exposants % 2 pour avoir que des 0 ou 1 
                # .count pour compter le nombre de facteurs pour ainsi avoir le nombre total d'exposant
        print("Liste des exposants : "+str(liste_exposant))
        print("")
        
        # Recherches de carré
        count+=1        
        if(1 not in liste_exposant): 
            carre=True
            return carre, x,count # On a trouvé 1 carré on peut dire conclure vers la factorisation
        #else: 
        #    pass
        M.append(liste_exposant)
        count+=1  
    print("")
    print("Matrice M : "+str(M))
    
    np.array(M)
    return(carre, np.transpose(M),count)


# Elimination de Gauss Algortihme  
def Elimination_Gauss(M):    
    count=0
    
    # Reconversion vers une liste, pour utiliser les fonctions index
    M=M.tolist()
    # Vecteur qui vérifie si il y a bien un pivot sur chaque ligne
    pivot_verif = [False]*len(M[0])
    count+=1
    
    for i in range(len(M)): 
        #Définissions des lignes de la matrice
        ligne = M[i]
        count+=1
  
        # Recherche d'un pivot en parcourant la ligne 
        #for num in ligne:
        for m_i in range(0,len(ligne)):
            # On a trouvé un pivot
            count+=2
            if(ligne[m_i]==1):
                # On garde l'indice de la ligne du pivot, et on indique sur notre vecteur verifie True
                indice=ligne.index(ligne[m_i])
                count+=1
                pivot_verif[indice] = True
                count+=1
                
                for j in chain(range(0,i),range(i+1,len(M))): 
                # Recherche d'autres 1 sur les même colonnes 
                    count+=2
                    if(M[j][indice] == 1):
                        
                        for i in range(len(M[j])):
                            count+=1
                            M[j][i]=(M[j][i] + ligne[i])%2
                            count+=3
                break
    
    M = np.transpose(M)
    count+=1
        
    eq_param = [] # représente la liste des équations qui seront dépendants dans M qui va permettre 
    # de trouver le vect 
    
    # Recherche de colonnes libres qui vont devenir des lignes
    for i in range(len(pivot_verif)): 
        count+=2
        if(pivot_verif[i]== False):
            ligne_libre = [M[i],i]
            count+=1
            eq_param.append(ligne_libre) # les lignes ou on va trouver des solutions
            count+=1
    
    # not eq_param
    count+=1
    if(eq_param==[]):
        return("Pas de Solutions ! Augmenter le nombre de nombres friables !")
        
    print("Solutions Trouvés : {}".format(len(eq_param)))
    
    # retourne le système d'équations à résoudre pour le vect, vecteur pivot verifie, matrice M
    return eq_param,pivot_verif,M,count


# Algorithme Vect va trouver les vects dans les équations paramétrés, (équations qui n'ont pas de pivots) 
def Vect(eq_param,M,pivot_verif,K=0):
    count=0
    vect=[] 
    indices = []
    ligne_libre = eq_param[K][0] # ligne multiple d'une autre ligne, L1 =2*L2 par exemple
    count+=1
    for i in range(0,len(ligne_libre)):
        count+=2
        if(ligne_libre[i] == 1): 
            indices.append(i)
            count+=1
            
    # ligne avec les 1 dans les mêmes colonnes = dépendante / liée   
    for ligne in range(0,len(M)): 
        for index in indices:
            count+=3
            if(M[ligne][index] == 1 and pivot_verif[ligne]):
                
                vect.append(ligne)
                count+=1
                break
            
    vect.append(eq_param[K][1])  
    count+=1     
    return vect,count
#%%
##################################################################
## -------------------  FACTORISATION    -------------------------
##################################################################
def Factorisation(vect,liste_num_friables,a_liste,N):
    count=0
    nb_solution=[]
    for i in vect:
        nb_solution.append(liste_num_friables[i])
        count+=1
    print("nb solution y : "+str(nb_solution))

    nb_value=[]
    for i in vect:
        nb_value.append(a_liste[i])
        count+=1
    print("nb value x : "+str(nb_value))
    
    # On multiplie tous les facteurs premiers
    ysol = 1
    for n in nb_solution:
        ysol =ysol*n
        
    print("y sol : "+str(ysol))
        
    # On multiplie tous les x    
    x = 1
    for n in nb_value:
        x=x*n
    
    print("x sol : "+str(x))
    # Racine Carré de ysol
    y,c = large_sqrt(ysol)
    count+=c
    print("large sqrt y : "+str(y))
    
    facteur = np.gcd(x-y,N)
    return facteur,count
#%%
liste_facteur=[]
count_Legendre=[]    
    

# Crible Quadratique, Version Polynome Simple, I l'intervalle du crible
def Crible_Quadratique(n,B,Intervalle,count):
    
    # facteur de tolérance, Nombre de nombre premiers dans la base de facteurs +1 = nombres d'équatinos
    global T
    T=1


    ###############################################
    #-------- DEBUT DU CRIBLE QUADRATIQUE  --------
    ###############################################
    #print("")
    #print("----------------------")
    #print("- Crible Quadratique -")
    #print("----------------------")
    #print("")
    # Test si le nombre est premier on return N
    count+=1+NombrePremier(n)[1]
    if(NombrePremier(n)[0]):
        liste_facteur.append(n)
        return liste_facteur,count 
    
    # Si la racine carré de N est un entier,peut être erroné,on retourne la racine carré pr un très grand nombre
    count+=1
    
    
    if(len(str(n))<20):
        float_verif=str(np.sqrt(n))
        # On vérifie que la longueur max est le chiffre juste apres la virgule qui est 0, par ex: a=str(2.0=), a[2]=0 et len(a)=3
        index_virgule=float_verif.index('.')
        if(float_verif[index_virgule + 1 ]=='0'):
            if(index_virgule+1)==len(float_verif)-1:        
                #if(type(np.sqrt(n))==int):
                    print("Racine Trouvé : "+str(int(np.sqrt(n))))
                    if(NombrePremier(np.sqrt(n))[0]==True):
                        count+=large_sqrt(n)[1]
                        liste_facteur.append(large_sqrt(n)[0])            
                        liste_facteur.append(large_sqrt(n)[0])
                        count+=2
                        return liste_facteur,count
                    else: 
                        return Crible_Quadratique(int(np.sqrt(n)),B,Intervalle,count)
                        
    
    

    print("")
    print("------------------")    
    print("1. INITIALISATION")
    print("------------------")
    print("")
    
    ###############################################
    #-------------- INITIALISATION  ---------------
    ###############################################
    #B=Seuil_Friabilite(n)
    count+=1
    
    print("Seuil de Friabilité B : {}".format(B))
    base_facteur = Base_Facteurs(n,B)[0]#generates a B-smooth factor base
    count+= Base_Facteurs(n,B)[1]
    print("Génération de la base de facteurs premiers : {}".format(base_facteur))

    # Nombre de nombres premiers dans la base de facteurs 
    Tho = len(base_facteur)
    
    print("On recherche  {} relations qui sont {}-friables".format(Tho+T,B))
    
    print("")
    print("------------")    
    print("2. CRIBLAGE")
    print("------------")
    print("")
    
    ###############################################
    #---------------- CRIBLAGE  -------------------
    ###############################################
    
    liste_num_friables,a_liste,index,c = Criblage(base_facteur, n,Intervalle)
    count+=c
    # On utilise l'algorithme de Tonelli Shanks pour trouver des racines/solutions à r²=a%p
    # Ainsi on trouvera toutes les relations qui sont B-Friables
    
    print("On a trouvé {} nombres qui sont B-friables :".format(len(liste_num_friables)))
   
    print("Liste des nombres friables : "+str(liste_num_friables))
    print("")
    print("Racine de n : "+str(large_sqrt(n)[0]))
    print("a liste : "+str(a_liste))
    print("")
    print("indice : "+str(index))
    # Si on a pas assez de nombres friables, on doit augmenter l'intervalle du criblage ou 
    # la le nombre de nombres premiers dans la base de facteurs
    count+=1
    if(len(base_facteur)>len(liste_num_friables)):
        return("Pas assez de nombres friables ! Il faut augmenter l'Intervalle du criblage ou La taille de base facteur !") 
    
    
    print("")
    print("--------------------")    
    print("3. ALGEBRE LINEAIRE")
    print("--------------------")
    print("")
    
    ###############################################
    #------------ ALGEBRE LINEAIRE  ---------------
    ###############################################
    
    # Génère la matrice M % 2 à partir des relations
    carre, M_transpose,c = Matrice(liste_num_friables,base_facteur)
    count+=c 

    print("carre : "+str(carre))
    print("Matrice Transpose : "+str(M_transpose))
    

    count+=1
    if(carre == True):
        x = liste_num_friables.index(M_transpose)
        count+=1
        facteur = np.gcd(a_liste[x]+int(np.sqrt(M_transpose)),n)
        count+=1
        print("On a trouvé un carré!")
        count+=1+NombrePremier(int(facteur))[1]
        
        
        if(NombrePremier(int(facteur))[0]==True):
                liste_facteur.append(int(facteur))
                count+=1
                
                # tous les facteurs sont premiers on a finis
                count+=1+NombrePremier(int(n//facteur))[1]
                if(NombrePremier(int(n//facteur))[0]==True):
                    liste_facteur.append(int(n//facteur))
                    count+=1
                    
                    print("Listes des facteurs non trivials :")
                    return liste_facteur,count
                
                # il reste encore à décomposer en nombre premier
                else:
            
                    return Crible_Quadratique(int(n//facteur),B,Intervalle,count)
                               
        # Le facteur trouvé n'est pas premier
        else:
                autre_facteur=int(n//facteur)
                count+=1+NombrePremier(autre_facteur)[1]
                if(NombrePremier(autre_facteur)[0]==True):
                    liste_facteur.append(autre_facteur)
                    count+=1
                    return Crible_Quadratique(int(n//autre_facteur),B,Intervalle,count)
                
                else:
                    return Crible_Quadratique(facteur,B,Intervalle,count)
    
    
    print("Elimination de Gauss")
    eq_param, pivot_verif, M, c = Elimination_Gauss(M_transpose) 
    #Résout matrice pour MX= 0, trouve le carré parfait
    count+=c
    print("Equations Paramétrées : "+str(eq_param))
    print("Pivot Verifie : "+str(pivot_verif))
    
    vect,c = Vect(eq_param,M,pivot_verif,0)
    print("Vect : "+str(vect))
    count+=c
    #print("")
 
    ###############################################
    #------------- FACTORISATION   ---------------#
    ###############################################
    print("")
    print("-----------------")    
    print("4. FACTORISATION")
    print("-----------------")
    print("")

    #Résolution des congruences au carré pour avoir les facteurs non trivials
    facteur,c = Factorisation(vect,liste_num_friables,a_liste,n)
    count+=c
    
    for i in range(1,len(eq_param)):
        count+=3
        if (facteur == 1 or facteur == n):
            
            print("Faux ! Essayez d'autres solutions vect!")
            vect,c = Vect(eq_param,M,pivot_verif,i)
            count+=c
            facteur,c = Factorisation(vect,liste_num_friables,a_liste,n)
            count+=c
        else:
            count+=1+NombrePremier(int(facteur))[1]
            if(NombrePremier(int(facteur))[0]==True):
                liste_facteur.append(int(facteur))
                count+=1
                
                # tous les facteurs sont premiers on a finis
                count+=1+NombrePremier(int(n//facteur))[1]
                if(NombrePremier(int(n//facteur))[0]==True):
                    liste_facteur.append(int(n//facteur))
                    count+=1
                    
                    #print("Listes des facteurs non trivials :")
                    print("")
                    print("-------------------------------")
                    print("")
                    return liste_facteur,count
                
                # il reste encore à décomposer en nombre premier
                else:
            
                    return Crible_Quadratique(int(n//facteur),B,Intervalle,count)
                               
            # Le facteur trouvé n'est pas premier
            else:
                autre_facteur=int(n//facteur)
                count+=1+NombrePremier(autre_facteur)[1]
                if(NombrePremier(autre_facteur)[0]==True):
                    liste_facteur.append(autre_facteur)
                    count+=1
                    return Crible_Quadratique(int(n//autre_facteur),B,Intervalle,count)
                
                else:
                    return Crible_Quadratique(facteur,B,Intervalle,count)
      
    return(liste_facteur,count,"Aucun Facteurs Non Triviaux !")
#%%
def CribleQuadratique_Chrono(n,B,I):
    time1=time.time() 
    resultat=Crible_Quadratique(n,B,I,0)
    time2=time.time() 
    time3=time2-time1
    return resultat,time3

#fact,count=Crible_Quadratique(1811706971,1000,0)
#print("Listes des facteurs : "+str(fact))
#print("compteur: "+str(len(count_Legendre)+count))

#n=146182562237
#n=1811706971
n=264839967043414254127


#n=16
#B=10
#I=10
#print(large_sqrt(1811706971*7841* 4421*4217))

#n=4
#n=16

n=264839967043414254127

Intervalle = 500000
B=3000
#n=10153331*373587883*76695841 #290919706205487829572593
#Intervalle=1000000
#B=4500
n=44899460164117431289
Intervalle=390000
B=1000

#print(Crible_Quadratique(146182562237,20000,0))
#print(Crible_Quadratique(1811706971,1000,0))
#print(Crible_Quadratique(n,B,Intervalle,0))



n=703253 
B=90
Intervalle=500


n=149*101*137
B=50
Intervalle=500

n=87463
B=19
Intervalle=120



print("")
print("----------------------")
print("- Crible Quadratique -")
print("----------------------")
print("")
print("Nombre à Factoriser n = "+str(n))
print("")
result=CribleQuadratique_Chrono(n,B,Intervalle)
print("Liste des facteurs : "+str(result[0][0]))
print("Nombres d'Opérations : "+str(result[0][1]))
print("Temps : "+str(result[1]))

#%%
"""
Bq=30
Intervalle=1000
def plot_cb(start,end,B,I):
    
    yop_p1=[]
    ychrono_p1=[]
    
    for i in range(start,end,10000):
        print("Iteration : "+str(i))
        result=CribleQuadratique_Chrono(i,B,I)
        print(result)
        
        yop_p1.append(result[0][1])
        ychrono_p1.append(result[1])
        print("")
        
    return [yop_p1,ychrono_p1]


import matplotlib as plt
import matplotlib.pyplot as plt


start=5
end=1000005

res_cb=plot_cb(start,end,Bq,Intervalle)
yop_cb=res_cb[0]
ychrono_cb=res_cb[1]
      
x= [i for i in range(start,end,10000)]


plt.plot(x, yop_cb,label = "Crible Quadratique")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Opérations nécessaires')
axes.set_title('Crible Quadratique')
plt.show()

plt.plot(x, ychrono_cb,label = "Crible Quadratique")
axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de '+str(start)+' à '+str(end))
axes.set_ylabel('Temps en secondes')
axes.set_title('Crible Quadratique')
plt.show()
"""

