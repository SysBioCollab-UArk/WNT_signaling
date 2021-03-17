from pysb import *
import numpy as np
import matplotlib.pyplot as plt

Model()

#monomers
Monomer("Bcat",["tcf4","dcplx","bcat","ub"],{"ub":["U","D"]})#Bcat
Monomer("Gli2")#Gli2 protein
Monomer("gGli2")#Gli2 gene
Monomer("mGli2")#Gli2 mRNA
Monomer("Wnt")#Wnt
Monomer("Btrcp")#Btrcp
Monomer("Tcf4")#Tcf4
Monomer("Smad3")#Smad3
Monomer("Dcplx")#Dcplx
Monomer("Rec")#Rec
Monomer("Dvl")#Dvl
Monomer("gPthlh")#Pthlh gene
Monomer("mPthlh")#Pthlh mRNA
Monomer("Pthrp")#Pthrp protein

print(model)
print(model.monomers)
