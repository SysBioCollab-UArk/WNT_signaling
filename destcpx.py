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
Monomer('Apc', ['axin', 'aa15', 'aa20', 'state'], {'state': ['u', 'p']})
Monomer('Bcat', ['top', 'bottom', 'nterm', 'btrcp', 'state'], {'state': ['x', 'ub'], 'nterm' : ['u', 'p1', 'p2']}) # p1 by ck1a, p2 by gsk3
Monomer('Btrcp', ['bcat'])

# Initials
Parameter('Axin_0', 10)
Parameter('Gsk3_0', 50)
Parameter('Ck1a_0', 50)
Parameter('Apc_0', 50)
Parameter('Bcat_0', 100)
Parameter('Btrcp_0', 50)

Initial(Axin(bcat=None, gsk3=None, ck1a=None, apc=None), Axin_0)
Initial(Gsk3(axin=None), Gsk3_0)
Initial(Ck1a(axin=None), Ck1a_0)
Initial(Apc(axin=None, aa15=None, aa20=None, state='u'), Apc_0)
Initial(Bcat(top=None, bottom=None, nterm='u', btrcp=None, state='x'), Bcat_0)
Initial(Btrcp(bcat=None), Btrcp_0)

#observables
Observable('Axin_tot', Axin())   # total amount of axin
Observable('Gsk3_tot', Gsk3())   # total amount of gsk3
Observable('Ck1a_tot', Ck1a())   # total amount of ck1a
Observable('Apc_aa20u', Apc(state='u'))
Observable('Apc_aa20p', Apc(state='p'))
Observable('Bcat_u', Bcat(nterm='u'))
Observable('Bcat_p1', Bcat(nterm='p1'))
Observable('Bcat_p2', Bcat(nterm='p2'))
Observable('Bcat_ub', Bcat(state='ub'))
Observable('Btrcp_tot', Btrcp())

# Rate constants
Parameter('kf_axin_ck1a',1)
Parameter('kf_axin_gsk3',1)
Parameter('kf_axin_apc',1)
Parameter('kf_bcat_dtcpx', 1)
Parameter('kf_bcat_apc', 1)
Parameter('kf_bcat_phos_gsk3', 1)
Parameter('kf_bcat_phos_ck1a', 1)
Parameter('kf_phos_ck1a', 2)
Parameter('kf_bcat_binds_apc', 1)
Parameter('k_btrcp_binds_bcat', 1)
Parameter('k_bcat_ubiq', 1)
Parameter('k_bcat_release', 2)
# Rules

# Axin binding rules

Rule('axin_binds_ck1a', Axin(ck1a=None) + Ck1a(axin=None) >> Axin(ck1a=1) % Ck1a(axin=1), kf_axin_ck1a)
Rule('axin_binds_gsk3', Axin(gsk3=None) + Gsk3(axin=None) >> Axin(gsk3=1) % Gsk3(axin=1), kf_axin_gsk3)
Rule('axin_binds_apc', Axin(apc=None) + Apc(axin=None) >> Axin(apc=1) % Apc(axin=1), kf_axin_apc)

# Bcat into dest comp

Rule('Bcat_binds_dtcpx', Bcat(top=None) + Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=ANY) >> \
     Bcat(top=1) % Axin(bcat=1, ck1a=ANY, gsk3=ANY, apc=ANY), kf_bcat_dtcpx )
Rule('Bcat_binds_aa15', Bcat(top=ANY,bottom=None) % Apc(aa15=None) >> Bcat(top=ANY, bottom=1) % Apc(aa15=1), kf_bcat_apc)

# is apc % bcat at aa15 necessary? not sure when it unbinds but has to unbind


# Bcat phosphorilated by gsk3 and ck1a

# Rule('bcat_p_ck1a', Bcat(nterm='u', top=ANY) >> Bcat(nterm='p1', top=ANY), kf_bcat_phos_ck1a) # ck1a phos first
# Rule('bcat_p_gsk3', Bcat(nterm='p1', top=ANY) >> Bcat(nterm='p2', top=ANY), kf_bcat_phos_gsk3) # then gsk3b phos after

# APC phosphorilated by ck1a
#should this be reversible or onedirectional?---AF

# Rule('apc_p_ck1a', Apc(state='u', axin=ANY) >> Apc(state='p', axin=ANY), kf_phos_ck1a)

# phosphoralation forces bcat detach from axin

# Rule('bcat_binds_apc', Bcat(top=1, nterm='p2') % Axin(bcat=1, ck1a=ANY, gsk3=ANY, apc=2) % \
#      Apc(aa20=None, axin=2, state='p') >> Bcat(top=1, nterm='p2') % Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=2) \
#      % Apc(aa20=1, axin=2, state='p'), kf_bcat_binds_apc)

# apc is still bound to axin and bcat and apc are both phos


#  Btrcp binds bcat

# Rule('btrcp_binds_bcat', Bcat(top=1, nterm='p2', btrcp=None) % Apc(state='p', aa20=1, axin=ANY) + Btrcp(bcat=None) >> \
#      Bcat(top=1, nterm='p2', btrcp=2) % Apc(state='p', aa20=1, axin=ANY) % Btrcp(bcat=2), k_btrcp_binds_bcat)
#
# # should axin=any be included?
#
#
# #  Btrcp ubiquitinates bcat
#
# Rule('Bcat_ubiq', Bcat(btrcp=1, state='x', top=ANY) % Btrcp(bcat=1) >> \
#      Bcat(btrcp=1, state='ub', top=ANY) % Btrcp(bcat=1), k_bcat_ubiq)
#
# # Apc release bcat by dephos
#
# Rule('apc_release_bcat', Apc(state='p', aa20=1, axin=ANY) % Bcat(top=1, state='ub', btrcp=ANY) >> \
#      Apc(state='u', aa20=None, axin=ANY) + Bcat(top=None, state='ub', btrcp=ANY), k_bcat_release)

#running simulations
# tspan = np.linspace(0, 1, 101)
# sim = ScipyOdeSimulator(model, tspan, verbose=True)
print(model.monomers)
print()
print(model.rules)
print()
generate_equations(model,verbose=True)
for i,sp in enumerate(model.species):
     print(i,sp)
print()
for i,rxn in enumerate(model.reactions):
     print(i,rxn)
