from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from pysb.bng import generate_equations

Model()

# monomer
Monomer('Axin', ['bcat', 'gsk3', 'ck1a', 'apc'])
Monomer('Gsk3', ['axin'])
Monomer('Ck1a', ['axin'])
Monomer('Apc', ['axin', 'aa15', 'aa20', 'state'], {'state': ['u', 'p1', 'p2']})
Monomer('Bcat', ['top', 'bottom', 'nterm', 'state'], {'state': ['x', 'ub'], 'nterm': ['u', 'p1', 'p2']}) # p1 by ck1a, p2 by gsk3
Monomer('Btrcp', ['bcat'])

# Initials
Parameter('Axin_0', 100)
Parameter('Gsk3_0', 50)
Parameter('Ck1a_0', 50)
Parameter('Apc_0', 50)
Parameter('Bcat_0', 100)
Parameter('Btrcp_0', 50)

Initial(Axin(bcat=None, gsk3=None, ck1a=None, apc=None), Axin_0)
Initial(Gsk3(axin=None), Gsk3_0)
Initial(Ck1a(axin=None), Ck1a_0)
Initial(Apc(axin=None, aa15=None, aa20=None, state='u'), Apc_0)
Initial(Bcat(top=None, bottom=None, nterm='u', state='x'), Bcat_0)
Initial(Btrcp(bcat=None), Btrcp_0)

#observables
# Observable('Axin_tot', Axin())   # total amount of axin
# Observable('Gsk3_tot', Gsk3())   # total amount of gsk3
# Observable('Ck1a_tot', Ck1a())   # total amount of ck1a
# Observable('Apc_aa20u', Apc(state='u'))
# Observable('Apc_aa20p', Apc(state='p'))
# Observable('Bcat_u', Bcat(nterm='u'))
# Observable('Bcat_p1', Bcat(nterm='p1'))
# Observable('Bcat_p2', Bcat(nterm='p2'))
# Observable('Bcat_ub', Bcat(state='ub'))
# Observable('Btrcp_tot', Btrcp())
Observable('ck1a_axin', Axin(ck1a=ANY))
Observable('gsk3_axin', Axin(gsk3=ANY))
Observable('apc_axin', Axin(apc=ANY))
# Observable('axin_free', Axin(bcat=None, gsk3=None, ck1a=None, apc=None))
# Observable('bcat_free', Bcat(top=None, bottom=None))
# Observable('bcat_axin', Bcat(top=ANY, bottom=None))
# Observable('bcat_axin_apc', Bcat(top=ANY, bottom=ANY))
Observable('dcplx', Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=ANY))
Observable('bcat_dcplx', Axin(bcat=ANY, ck1a=ANY, gsk3=ANY, apc=ANY))
Observable('bcat_p1', Bcat(top=ANY, nterm='p1'))
Observable('bcat_p2', Bcat(top=ANY, nterm='p2'))
Observable('apc_p1', Apc(state='p1'))
Observable('apc_p2', Apc(state='p2'))
Observable('bcat_apc_aa15', Apc(aa15=ANY, aa20=None))
Observable('bcat_apc_aa20', Apc(aa15=None, aa20=ANY))
Observable('bcat_apc_both', Apc(aa15=ANY, aa20=ANY))
Observable('bcat_ub_apc', Bcat(top=1, bottom=None, state='ub') % Apc(aa20=1))
Observable('bcat_ub_btrcp', Bcat(state='ub') % Btrcp())
Observable('bcat_ub_total', Bcat(state='ub'))

OBS = [
    ['ck1a_axin', 'gsk3_axin', 'apc_axin'],
    ['dcplx', 'bcat_dcplx'],
    ['bcat_p1', 'bcat_p2'],
    ['apc_p1', 'apc_p2'],
    ['bcat_apc_aa15', 'bcat_apc_aa20', 'bcat_apc_both'],
    ['bcat_ub_apc', 'bcat_ub_btrcp', 'bcat_ub_total']
]

# Rate constants
Parameter('kf_axin_ck1a', 1)
Parameter('kr_axin_ck1a', 10)
Parameter('kf_axin_gsk3', 1)
Parameter('kr_axin_gsk3', 100)
Parameter('kf_axin_apc', 1)
Parameter('kr_axin_apc', 50)
Parameter('kf_bcat_dtcpx', 100)
Parameter('kr_bcat_dtcpx', 0.1)
Parameter('kf_bcat_apc', 100)
Parameter('kr_bcat_apc', 0.1)
Parameter('kf_bcat_phos_gsk3', 100)
# Parameter('kr_bcat_phos_gsk3', 1)
Parameter('kf_bcat_phos_ck1a', 10)
# Parameter('kr_bcat_phos_ck1a', 1)
Parameter('kf_apc_phos_ck1a', 20)
# Parameter('kr_apc_phos_ck1a', 2)
Parameter('kf_apc_phos_gsk3', 200)
Parameter('k_dephos', 0.1)
Parameter('kf_bcat_binds_apc', 100)
Parameter('k_btrcp_binds_bcat', 1)
Parameter('k_bcat_ubiq', 1)
Parameter('k_bcat_release', 0.1)
Parameter('k_bcat_deg', 1)
Parameter('k_apc_dephos', 10)

# Rules

# Axin binding rules

Rule('axin_binds_ck1a', Axin(bcat=None, ck1a=None) + Ck1a(axin=None) | Axin(bcat=None, ck1a=1) % Ck1a(axin=1),
     kf_axin_ck1a, kr_axin_ck1a)
Rule('axin_binds_gsk3', Axin(bcat=None, gsk3=None) + Gsk3(axin=None) | Axin(bcat=None, gsk3=1) % Gsk3(axin=1),
     kf_axin_gsk3, kr_axin_gsk3)
Rule('axin_binds_apc', Axin(bcat=None, apc=None) + Apc(axin=None) | Axin(bcat=None, apc=1) % Apc(axin=1),
     kf_axin_apc, kr_axin_apc)

# Bcat into dest comp

Rule('Bcat_binds_dtcpx',
     Bcat(top=None, bottom=None, nterm='u') +
     Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=ANY) % Apc(aa20=None) |
     Bcat(top=1, bottom=None, nterm='u') % Axin(bcat=1, ck1a=ANY, gsk3=ANY, apc=ANY) % Apc(aa20=None),
     kf_bcat_dtcpx, kr_bcat_dtcpx)

# We think Bcat can only be phosphorylated and dephosphorylated when bound to Axin, not to APC
Rule('bcat_p_ck1a', Bcat(top=1, nterm='u') % Axin(bcat=1) | Bcat(top=1, nterm='p1') % Axin(bcat=1),
     kf_bcat_phos_ck1a, k_dephos)  # ck1a phos first
# Rule('bcat_unp1',   Bcat(nterm='p1') >> Bcat(nterm='u'), k_dephos)
Rule('bcat_p_gsk3', Bcat(top=1, nterm='p1') % Axin(bcat=1) | Bcat(top=1, nterm='p2') % Axin(bcat=1),
     kf_bcat_phos_gsk3, k_dephos)  # then gsk3b phos after
# Rule('bcat_unp2',   Bcat(nterm='p2') >> Bcat(nterm='p1'), k_dephos)

Rule('Bcat_binds_aa15', Axin(bcat=1) % Bcat(top=1, bottom=None) % Apc(state='u', aa15=None) |
     Axin(bcat=1) % Bcat(top=1, bottom=2) % Apc(state='u', aa15=2), kf_bcat_apc, kr_bcat_apc)

# is apc % bcat at aa15 necessary? not sure when it unbinds but has to unbind
# It aslso binds to the aa20 but when it is phosphorilated .According to the paper
# As discussed below, we have suggested that the reason APC contains two different types of b-catenin binding motifs is because they may function differently within the destruction
# complex.

# Bcat phosphorylated by gsk3 and ck1a

# According to a paper in the absence of WNT signals the destruction is activated and a balanced between unphos and phos B-catenin

# APC phosphorylated by ck1a
# should this be reversible or onedirectional?---AF

Rule('apc_p_ck1a', Apc(aa15=1, state='u', axin=ANY) % Ck1a() % Gsk3() % Bcat(bottom=1, nterm='p2') >>
     Apc(aa15=1, state='p1', axin=ANY) % Ck1a() % Gsk3() % Bcat(bottom=1, nterm='p2'), kf_apc_phos_ck1a)
Rule('apc_unp1', Apc(state='p1', aa20=None) >> Apc(state='u', aa20=None), k_dephos)

Rule('apc_p_gsk3', Apc(aa15=1, state='p1', axin=ANY) % Ck1a() % Gsk3() % Bcat(bottom=1, nterm='p2') >>
     Apc(aa15=1, state='p2', axin=ANY) % Ck1a() % Gsk3() % Bcat(bottom=1, nterm='p2'), kf_apc_phos_gsk3)
Rule('apc_unp2', Apc(state='p2', aa20=None) >> Apc(state='p1', aa20=None), k_dephos)

# phosphorylation forces bcat detach from axin

Rule('bcat_binds_apc', Bcat(top=1, nterm='p2', bottom=3) % Axin(bcat=1) % Ck1a() % Gsk3() %
     Apc(aa20=None, state='p2', aa15=3) >> Bcat(top=2, nterm='p2', bottom=None) % Axin(bcat=None) % Ck1a() % Gsk3() %
     Apc(aa20=2, state='p2', aa15=None), kf_bcat_binds_apc)

# apc is still bound to axin and bcat and apc are both phos

#  Btrcp binds bcat

Rule('btrcp_binds_bcat',
     Bcat(top=1, nterm='p2', bottom=None) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) + Btrcp(bcat=None) >>
     Bcat(top=1, nterm='p2', bottom=2) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) % Btrcp(bcat=2),
     k_btrcp_binds_bcat)

# should axin=any be included?

#  Btrcp ubiquitinates bcat
Rule('Bcat_ubiq',
     Bcat(top=1, bottom=2, state='x') % Apc(aa20=1) % Btrcp(bcat=2) >>
     Bcat(top=1, bottom=2, state='ub') % Apc(aa20=1) % Btrcp(bcat=2), k_bcat_ubiq)

# Apc release bcat by dephos
Rule('apc_release_bcat',
     Bcat(top=1, bottom=2, state='ub') % Apc(aa20=1) % Btrcp(bcat=2) >>
     Bcat(top=None, bottom=2, state='ub') % Btrcp(bcat=2) + Apc(aa20=None), k_bcat_release)

# APC degraded by proteosome
Rule('bcat_degradation', Bcat(top=None, bottom=2, state='ub') % Btrcp(bcat=2) >> Btrcp(bcat=None), k_bcat_deg)

# dephosphorylation of APC after Bcat release
Rule('apc_dephos', Apc(state='p2', aa20=None) >> Apc(state='u', aa20=None), k_apc_dephos)

#running simulations
tspan = np.linspace(0, 1, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()
# for obs in model.observables:
     # if obs.name == 'bcat_free':
     #      plt.legend(loc=0)
     #      plt.figure()
for group in OBS:
    plt.figure()
    for obs_name in group:
        plt.plot(tspan, result.observables[obs_name], lw=2, label=obs_name)
    plt.legend(loc=0)
    plt.xlabel('time')
    plt.ylabel('concentration')

# print(len(model.species))
# for sp in model.species:
#      print(sp)

plt.show()

# print(model.monomers)
# print()
# print(model.rules)
# print()
# generate_equations(model,verbose=True)
# for i,sp in enumerate(model.species):
#      print(i,sp)
# print()
# for i,rxn in enumerate(model.reactions):
#      print(i,rxn)
