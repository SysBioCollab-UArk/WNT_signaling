from pysb import *
import numpy as np
import matplotlib.pyplot as plt

Model()

# monomers
Monomer('Bcat', ['tcf4', 'dcplx', 'bcat', 'ub', 'dvl'], {'ub': ['U', 'D']})  # Bcat
Monomer('Gli2', ['btrcp', 'gPthlh'])  # Gli2 protein
Monomer('gGli2', ['tcf4', 'smad3'])  # Gli2 gene
Monomer('mGli2')  # Gli2 mRNA
Monomer('Wnt', ['rec'])  # Wnt
Monomer('Btrcp', ['gli2', 'bcat'])  # Btrcp
Monomer('Tcf4', ['gGli2', 'gPthlh', 'bcat'])  # Tcf4
Monomer('Smad3', ['gGli2'])  # Smad3
Monomer('Dcplx', ['bcat'])  # Destruction complex
Monomer('Rec', ['wnt', 'dvl'])  # Receptor
Monomer('Dvl', ['rec', 'bcat'])  # Dvl
Monomer('gPthlh', ['gli2', 'tcf4'])  # Pthlh gene
Monomer('mPthlh')  # Pthlh mRNA
Monomer('Pthrp')  # Pthrp protein

# parameter

# rules

Rule('Bcat_u', Bcat(ub=D) + Dcplx(bcat=None) >> Bcat(ub=U) % Dcplx(bcat=1), kubiq)
Rule('Bcat_naked', Bcat(dcplx=2,dvl=None) + Dvl(bcat=None) >> Bcat(dcplx=None,dvl=3) % Dvl(bcat=3), knaked)
Rule('Bcat_free', Bcat(dcplx=None,dvl=3) % Dvl(bcat=3) >> Bcat(dcplx=None,dvl=None) + Dvl(bcat=None), kfree)
Rule('Bcat_tcf4', Bcat(dcplx=None,tcf4=None) + Tcf4(bcat=None) >> Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4), kpromg)
Rule('gGli2_prod', Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4,gGli2=5) + Smad3(gGli2=5) >> Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4,gGli2=5) + Smad3(gGli2=5) +mGli2, kgprod)
Rule('Dvl_bound', Wnt(rec=6) + Dvl(rec=None) | Wnt(rec=6) + Dvl(rec=7))