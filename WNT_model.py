from pysb import *
import numpy as np
import matplotlib.pyplot as plt

Model()

# monomers
Monomer('Bcat', ['tcf4', 'gsk3b', 'apc', 'state', 'btrcp'], {'state': ['x', 'ub']})  # Bcat
Monomer('Gli2', ['btrcp', 'gPthlh', 'state'], {'state': ['x', 'p', 'ub']})  # Gli2 protein
Monomer('gGli2', ['tcf4', 'smad3'])  # Gli2 gene
Monomer('mGli2')  # Gli2 mRNA
Monomer('Wnt', ['rec'])  # Wnt, rec is receptor complex - lrp5 and fzd
Monomer('Btrcp', ['b'])  # Btrcp
Monomer('Tcf4', ['gGli2', 'gPthlh', 'bcat'])  # Tcf4
Monomer('Smad3', ['gGli2'])  # Smad3
Monomer('Gsk3b', ['bcat', 'dvl'])  # represents gsk3, axin, ck1
Monomer('Rec', ['wnt', 'dvl'])  # Receptor
Monomer('Dvl', ['rec', 'gsk3b'])  # Dvl
Monomer('gPthlh', ['gli2', 'tcf4'])  # Pthlh gene
Monomer('mPthlh')  # Pthlh mRNA
Monomer('Pthrp')  # Pthrp protein
Monomer('Dkk1', ['rec'])   # wnt inhibitor at site
Monomer('Wif1', ['wnt'])   # wnt inhibitor at ligand
Monomer('Apc', ['bcat'])   # apc

# initials

Initial(Bcat(tcf4=None,apc=None,gsk3b=None,state='x',btrcp=None), Parameter('Bcat_0', 50))
Initial(Apc(bcat=None), Parameter('Apc_0', 500))
Initial(Gsk3b(bcat=None,dvl=None), Parameter('Gsk3b_0', 50))
Initial(Dvl(rec=1,gsk3b=None) % Rec(dvl=1,wnt=2) % Wnt(rec=2), Parameter('Dvl_rec_wnt_0', 50))
Initial(Gli2(btrcp=None,state='ub',gPthlh=None), Parameter('Gli2_0', 50))
Initial(Btrcp(b=None), Parameter('Btrcp_0', 50))






# parameter

Parameter('kubiq', 15)
Parameter('knaked', 2)
Parameter('kfree', 47)
Parameter('kpromg',12)
Parameter('kgprod', 14)
Parameter('kINH', .22)
Parameter('kapc', 53)
Parameter('kgsk', 70)

# rules

Rule('Bcat_u', Bcat(gsk3b=1,apc=2,state='x',btrcp=None) % Gsk3b(bcat=1) % Apc(bcat=2) + Gli2(btrcp=1,state='ub') % Btrcp(b=1) \
     >> Bcat(gsk3b=None,apc=None,state='ub',btrcp=1) % Btrcp(b=1) + Gsk3b(bcat=None) +Apc(bcat=None) + Gli2(btrcp=None,state='x'), \
kubiq)
Rule('Bcat_naked', Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=None) % Apc(bcat=2) + Dvl(gsk3b=None,rec=ANY) >> Bcat(gsk3b=1,apc=2) \
     % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Dvl(gsk3b=3,rec=ANY), knaked)
Rule('Bcat_apc', Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Dvl(gsk3b=3,rec=ANY) >> Bcat(gsk3b=1,apc=None) \
     % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) + Apc(bcat=None), kapc)
Rule('Bcat_gsk3b', Bcat(gsk3b=1,apc=None) % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) >> Bcat(gsk3b=None,apc=None) + Gsk3b(bcat=None,dvl=None) \
+ Dvl(gsk3b=None,rec=ANY), kgsk)
#Rule('Bcat_tcf4', Bcat(dcplx=None,tcf4=None) + Tcf4(bcat=None) >> Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4), kpromg)
#Rule('gGli2_prod', Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4,gGli2=5) + Smad3(gGli2=5) >> Bcat(dcplx=None,tcf4=4) % Tcf4(bcat=4,gGli2=5) + Smad3(gGli2=5) +mGli2, kgprod)
#Rule('Dvl_bound', Wnt(rec=6) + Dvl(rec=None) | Wnt(rec=6) + Dvl(rec=7))

from pysb.bng import generate_equations

generate_equations(model,verbose=True)