import numpy as np
import matplotlib as plt
import matplotlib.pyplot as plt

from AlgoNaif import Algo_Naif_Rec
from AlgoNaif import AlgoNaif_Chrono

from RhoPollard import Pollard_Rhos
from RhoPollard import Pollard_Rhos_Chrono

from P1Pollard import p_1_Pollard
from Lenstra import Lenstra
from CribleQuadratique import Crible_Quadratique


#%%
# Algo Naif

"""
nbPremier=2*3*5*7*11*2*2*5
nbPremier=146182562237
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


#%%
# Rho Pollard
"""
print("--------------------------------------------------------")
print("--              ALGORITHME RHO POLLARD                --")
print("--------------------------------------------------------")
print("")
n=146182562237
print("Test Algo RhoPollard")
print("Valeur n : "+str(n))
resultat_Rho=Pollard_Rhos(n,0)
print("Decomposition en Facteur Premier : "+str(resultat_Rho[0]))
print("Compteur : "+str(resultat_Rho[1]))
print("")

print("--------------------------------------------------------")
print("--            FIN ALGORITHME RHO POLLARD              --")
print("--------------------------------------------------------")
"""

#%%
#P-1 Pollard
"""
print("")
print("--------------------------------------------------------")
print("--              ALGORITHME P-1 POLLARD                --")
print("--------------------------------------------------------")
n=146182562237
print("List of Factor of "+str(n)+" is : "+str(p_1_Pollard(n))) 
print("")
print("--------------------------------------------------------")
print("--            FIN ALGORITHME P-1 POLLARD              --")
print("--------------------------------------------------------")
"""
#%%
# Lenstra
"""
print("")
print("--------------------------------------------------------")
print("--               ALGORITHME LENSTRA                   --")
print("--------------------------------------------------------")
print("")
n=2038452820458036
n=146182562237
#n=34665061
print("Facteur Non Trivial d : "+str(Lenstra(n)))
#4421
#7841
#4217

print("")
print("--------------------------------------------------------")
print("--             FIN ALGORITHME LENSTRA                 --")
print("--------------------------------------------------------")
print("")
"""

#%%
#Crible Quadratique
"""
print("")
print("--------------------------------------------------------")
print("--         DEBUT ALGORITHME CRIBLE QUADRATIQUE        --")
print("--------------------------------------------------------")
print("")

print(Crible_Quadratique(146182562237,20000))

print("")
print("--------------------------------------------------------")
print("---         FIN ALGORITHME CRIBLE QUADRATIQUE        ---")
print("--------------------------------------------------------")
print("")

"""

#%%
# Tracé du graphique de l'évolution du nombre d'opérations nécessaires en fonction de la méthode choisie
#n=146182562237


#start=146182562237
#n=146182562337

start=2
end=1000


"""
x= [i for i in range(start,end)]
y = [Algo_Naif_Rec(i,0)[1] for i in range(start,end)]

#y = [Pollard_Rhos(i,0)[1] for i in range(2,n)]

plt.plot(x, y,label = "Naif")

axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de 2 à n')
axes.set_ylabel('Opérations nécessaires')
axes.set_title('Graphique de l\'évolution du nombre d\'opérations nécessaires en fonction de la méthode choisie')
plt.show()



x= [i for i in range(start,end)]
y = [AlgoNaif_Chrono(i)[1] for i in range(start,end)]

axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de 2 à n')
axes.set_ylabel('Temps en secondes')
axes.set_title('Graphique de l\'évolution de la factorisation par la méthode Naïf en fonction du Temps en secondes')

plt.plot(x,y,label="Naif")
plt.show()

"""


#%%
#start=2
#end=1000

start=10000
end=100000
x= [i for i in range(start,end)]
y = [Pollard_Rhos(i,0)[1] for i in range(start,end)]

#y = [Pollard_Rhos(i,0)[1] for i in range(2,n)]

plt.plot(x, y,label = "Rho Pollard")

axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de 2 à n')
axes.set_ylabel('Opérations nécessaires')
axes.set_title('Graphique de l\'évolution du nombre d\'opérations nécessaires en fonction de la méthode choisie')
plt.show()



x= [i for i in range(start,end)]
y = [Pollard_Rhos_Chrono(i)[1] for i in range(start,end)]

axes = plt.gca()
axes.set_xlabel('Nombres à factoriser de 2 à n')
axes.set_ylabel('Temps en secondes')
axes.set_title('Graphique de l\'évolution de la factorisation par la méthode Rho Pollard en fonction du Temps en secondes')

plt.plot(x,y,label="Rho Pollard")
plt.show()










"""
for i in range(2,50):
    #print("i : "+str(i)+", Algo Naif : "+str(Algo_Naif_Rec(i,0)))
    print("i : "+str(i)+", Rho Pollard : "+str(Pollard_Rhos(i,0)))
"""
# Problème Algorithmique, la liste reste et ne se réinitialise pas 
    