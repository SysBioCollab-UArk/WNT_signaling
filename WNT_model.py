from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator


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

# initials

Initial(Bcat(tcf4=None,apc=None,gsk3b=None,state='x',btrcp=None), Parameter('Bcat_0', 50))
Initial(Apc(bcat=None), Parameter('Apc_0', 500))
Initial(Gsk3b(bcat=None,dvl=None), Parameter('Gsk3b_0', 50))
Initial(Dvl(rec=1,gsk3b=None) % Rec(dvl=1,wnt=2) % Wnt(rec=2), Parameter('Dvl_rec_wnt_0', 50))


# parameter

Parameter('kubiq', 15)
Parameter('knaked', 2)
Parameter('kfree', 47)
Parameter('kpromg',12)
Parameter('kgprod', 14)
Parameter('kINH', .22)

# rules

Rule('Bcat_u', Bcat(gsk3b=1,apc=2,state='x',btrcp=None) % Gsk3b(bcat=1) % Apc(bcat=2) + Gli2(btrcp=1,state='ub') \
     % Btrcp(b=1) >> Bcat(gsk3b=None,apc=None,state='ub',btrcp=1) % Btrcp(b=1) + Gsk3b(bcat=None) +Apc(bcat=None) \
     + Gli2(btrcp=None,state='x'), kubiq)
Rule('Bcat_naked', Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=None) % Apc(bcat=2) + Dvl(gsk3b=None,rec=ANY) >> \
     Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Dvl(gsk3b=3,rec=ANY), knaked)
Rule('Bcat_apc', Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Dvl(gsk3b=3,rec=ANY) >> \
     Bcat(gsk3b=1,apc=None) % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) + Apc(bcat=None), kapc)
Rule('Bcat_gsk3b', Bcat(gsk3b=1,apc=None) % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) >> Bcat(gsk3b=None,apc=None) \
     + Gsk3b(bcat=None,dvl=None) + Dvl(gsk3b=None,rec=ANY), kgsk)
Rule('Bcat_tcf4', Bcat(gsk3b=None,apc=None,tcf4=None) + Tcf4(bcat=None) >> Bcat(gsk3b=None,apc=None,tcf4=1) % \
     Tcf4(bcat=1), ktcf4)
Rule('gGli2', Bcat(gsk3b=None,apc=None,tcf4=1) % Tcf4(bcat=1,gGli2=None) + Smad3(gGli2=None) + gGli2(tcf4=None,smad3=None) \
     >> Bcat(gsk3b=None,apc=None,tcf4=1) % Tcf4(bcat=1,gGli2=2) % Smad3(gGli2=1) % gGli2(tcf4=2,smad3=1) + mGli2, kgGli2)
Rule('gPthlh', Bcat(gsk3b=None,apc=None,tcf4=1) % Tcf4(bcat=1,gPthlh=None) + Gli2(gPthlh=None,state='x') + gPthlh(gli2=None,tcf4=None) \
     >> Bcat(gsk3b=None,apc=None,tcf4=1) % Tcf4(bcat=1,gPthlh=2) % Gli2(gPthlh=1,state='x') % gPthlh(gli2=1,tcf4=2) + mPthlh(), kpthlh)
Rule('Pthrp', mPthlh() >> mPthlh() + Pthrp(), kpthrp)
Rule('gli2_p', mGli2() >> mGli2() + Gli2(state='x'), kgli2_p)
Rule('Rec_wnt', Wnt(rec=None) + Rec(wnt=None,dvl=None) + Dvl(rec=None) | Wnt(rec=1) % Rec(wnt=1,dvl=2) % Dvl(rec=2), krecw, krecwb)
# missing rules for competitive binding on receptor and wnt ligand, gli2 transforms


from pysb.bng import generate_equations
generate_equations(model,verbose=True)






# Create monomers of each of the proteins/compounds
Monomer('Wif1', ['wnt'])   # wnt inhibitor at ligand
Monomer('Apc', ['bcat'])   # apc

# Assign initial values to the starting monomers
Initial(Gli2(btrcp=None,state='ub',gPthlh=None), Parameter('Gli2_0', 50))
Initial(Btrcp(b=None), Parameter('Btrcp_0', 50))

# Create rules for chemical reactions
Rule('Rec_dkk1', Dkk1(rec=None) + Rec(wnt=None) | Dkk1(rec=1) % Rec(wnt=1), krecd, krecdb)
Rule('Wnt_inh', Wnt(rec=None) + Wif1(wnt=None) | Wnt(rec=1) % Wif1(wnt=1), kwnt, kwntb)

# Give each reaction rate a value
Parameter('kapc', 53)
Parameter('kgsk', 70)

# Observables
Observable('pthrp_tot', Pthrp())
Observable('gli2_free', Gli2(btrcp=None,gPthlh=None))

tspan=np.linspace(0,100,101)
sim=ScipyOdeSimulator(model,tspan,verbose=True)
traj=sim.run()
plt.plot(tspan,traj.observables['pthrp_tot'])

for obs in model.observables:
     plt.figure()
     plt.plot(tspan, traj.observables[obs.name], lw=2, label=obs.name)
     plt.legend(loc=0)
     plt.xlabel('Time')
     plt.ylabel('Concentration')