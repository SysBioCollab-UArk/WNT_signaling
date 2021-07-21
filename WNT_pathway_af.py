from pysb import *
import numpy as np
import matplotlib.pyplot as pl

Model()

#monomers
Monomer("Bcat", ["tcf4", "dcplx", "bcat", "ub"], {"ub": ["U", "D"]})#Bcat
Monomer("Gli2", ["btrcp", "ub", "P"], {"ub": ["U", "D"]}, {"P":["p+", "p-"]})#Gli2 protein it can also be phosphorylated"P"
Monomer("gGli2", ["smad3", "TCF4"])#Gli2 gene
Monomer("mGli2")#Gli2 mRNA
Monomer("Wnt", ["fzd", "LRP5"])#Wnt
Monomer("Btrcp", ["bcat"])#Btrcp But it only binds to B catenin when it is Ubiquitinated
Monomer("Tcf4", ["gGli2", "gPthlh", "Bcat"])#Tcf4
Monomer("Smad3", ["gGli2"])#Smad3
Monomer("Dcplx", ["bcat"])#Dcplx
Monomer("LRP5")#wnt coreceptor
Monomer("Dvl", ["DVl", "Dcplx", "LRP5"])#Dvl
Monomer("gPthlh")#Pthlh gene
Monomer("mPthlh")#Pthlh mRNA
Monomer("Pthrp")#Pthrp protein
Monomer("FZD", ["LRP5", "wnt","inh"])#frizzed receptor
Monomer("Inh", ["FZD"])# INHIBITORS include DKK1, WIF1,SFRP4
#Parameters
parameter('kubiq', 15)
parameter('knaked', 2)
parameter('kfree', 47)
parameter('kpromg',12)
parameter('kgprod', 14)
parameter('kINH'.22)

# rules

Rule('Bcat_u', Bcat(ub=D) + Dcplx(Bcat=None) >> Bcat(ub=U) % Dcplx(Bcat=1), kubiq)
Rule('Bcat_naked', Bcat(Dcplx=2, dvl=None) + Dvl(Dcpx=4) >> Bcat(Dcplx=None, dvl=3) % Dvl(bcat=3), knaked) #not really sure if DVL binds to Bcat so I CHANGED IT
Rule('Bcat_free', Bcat(dcplx=None,dvl=3) % Dvl(bcat=3) >> Bcat(dcplx=None,dvl=None) + Dvl(bcat=None), kfree)
Rule('Bcat_tcf4', Bcat(dcplx=None,tcf4=None) + Tcf4(bcat=None) >> Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4), kpromg)
Rule('gGli2_prod', Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4,gGli2=5) + Smad3(gGli2=5) >> Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4,gGli2=5) + Smad3(gGli2=5) +mGli2, kgprod)
Rule('Dvl_bound', Wnt(rec=6) + Dvl(rec=None) | Wnt(rec=6) + Dvl(rec=7))
#WNT INHIBITION
Rule('inh_FZD', FZD(inh=none, LRP5=3)% LRP5(FZD=5) >> FZD(inh=44, LRP5=5) % LRP5(FZD=5), kINH)

print(model)
print(model.monomers)
