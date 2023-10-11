from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from pysb.bng import generate_equations

# model based on van Kappel and Maurice 2017 (https://dx.doi.org/10.1111/bph.13922)
# rule verification was completed by LH and GB on 1/5/22

# 6 monomers + 30 rules (including the reverse directions)
# expands into 32 species + 63 reactions

# TODO: Modify rate constants so that more Bcat is degraded (very little is being degraded right now, 02/18/22)
# TODO: Should 'lithium' be a separate binding site on GSK3 or should lithium bind to the axin binding site?
Model()

# monomer
Monomer('Axin', ['bcat', 'gsk3', 'ck1a', 'apc'])
Monomer('Gsk3', ['axin', 'lithium'])
Monomer('Ck1a', ['axin'])
Monomer('Apc', ['axin', 'aa15', 'aa20', 'state'], {'state': ['u', 'p1', 'p2']})
Monomer('Bcat', ['top', 'bottom', 'nterm', 'state'], {'state': ['x', 'ub'], 'nterm': ['u', 'p1', 'p2']}) # p1 by ck1a, p2 by gsk3
Monomer('Btrcp', ['bcat'])
Monomer('Li', ['gsk3'])

# Initials
Parameter('Axin_0', 100)
Parameter('Gsk3_0', 50)
Parameter('Ck1a_0', 50)
Parameter('Apc_0', 50)
Parameter('Bcat_0', 100)
Parameter('Btrcp_0', 50)
Parameter('Li_0', 0)

Initial(Axin(bcat=None, gsk3=None, ck1a=None, apc=None), Axin_0)
Initial(Gsk3(axin=None, lithium=None), Gsk3_0)
Initial(Ck1a(axin=None), Ck1a_0)
Initial(Apc(axin=None, aa15=None, aa20=None, state='u'), Apc_0)
Initial(Bcat(top=None, bottom=None, nterm='u', state='x'), Bcat_0)
Initial(Btrcp(bcat=None), Btrcp_0)
Initial(Li(gsk3=None), Li_0)

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
Observable('bcat_p2', Bcat(top=ANY, nterm='p2'))  # We think phos of Bcat by GSK3 is "kinase activity" in Stambolic1996
Observable('apc_p1', Apc(state='p1'))
Observable('apc_p2', Apc(state='p2'))  # We think phos of Apc by GSK3 is "kinase activity" in Stambolic1996
Observable('bcat_apc_aa15', Apc(aa15=ANY, aa20=None))
Observable('bcat_apc_aa20', Apc(aa15=None, aa20=ANY))
Observable('bcat_apc_both', Apc(aa15=ANY, aa20=ANY))
Observable('bcat_ub_apc', Bcat(top=1, bottom=None, state='ub') % Apc(aa20=1))
Observable('bcat_ub_apc_btrcp', Bcat(top=1, bottom=2, state='ub') % Apc(aa20=1) % Btrcp(bcat=2))
Observable('bcat_ub_btrcp', Bcat(top=None, bottom=1, state='ub') % Btrcp(bcat=1))
Observable('bcat_ub_btrcp_total', Bcat(bottom=1, state='ub') % Btrcp(bcat=1))
Observable('bcat_ub_total', Bcat(state='ub'))
Observable('bcat_total', Bcat())
Observable('Li_total', Li())
Observable('Li_Gsk3', Li(gsk3=ANY))
###
Observable('GSK3_activity', Bcat(nterm='p2') + Apc(state='p2'))

OBS = [
    # ['ck1a_axin', 'gsk3_axin', 'apc_axin'],
    # ['dcplx', 'bcat_dcplx'],
    # ['bcat_p1', 'bcat_p2'],
    # ['apc_p1', 'apc_p2'],
    # ['bcat_apc_aa15', 'bcat_apc_aa20', 'bcat_apc_both'],
    # ['bcat_ub_apc', 'bcat_ub_apc_btrcp', 'bcat_ub_btrcp',
    #  'bcat_ub_btrcp_total', 'bcat_ub_total'],
    # ['bcat_ub_total', 'Li_total', 'Li_Gsk3'],
    ['bcat_total']
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
# Parameter('kr_bcat_phos_gsk3', 1)  # assuming (for now) that all dephos rates are equal to k_dephos
Parameter('kf_bcat_phos_ck1a', 10)
# Parameter('kr_bcat_phos_ck1a', 1)  # assuming (for now) that all dephos rates are equal to k_dephos
Parameter('kf_apc_phos_ck1a', 20)
# Parameter('kr_apc_phos_ck1a', 2)  # assuming (for now) that all dephos rates are equal to k_dephos
Parameter('kf_apc_phos_gsk3', 200)
Parameter('k_dephos', 0.1)
Parameter('kf_bcat_binds_apc', 100)
Parameter('kf_btrcp_binds_bcat', 0.01)
Parameter('kr_btrcp_binds_bcat', 1)
Parameter('k_bcat_ubiq', 0.1)#1
Parameter('k_bcat_release', 10)#1
Parameter('k_bcat_deg', 0.001)#0.1
Parameter('kf_gsk3_li', 1)
Parameter('kr_gsk3_li', 0.01)

# Rules

# STEP 1 ###########
# Axin binding rules
# We require beta-catenin to NOT be bound for these binding event to occur
# We assume beta-catenin can only bind when all three of Ck1a, Gsk3, and Apc are bound
Rule('axin_binds_ck1a', Axin(bcat=None, ck1a=None) + Ck1a(axin=None) |
     Axin(bcat=None, ck1a=1) % Ck1a(axin=1),
     kf_axin_ck1a, kr_axin_ck1a)
Rule('axin_binds_gsk3', Axin(bcat=None, gsk3=None) + Gsk3(axin=None) |
     Axin(bcat=None, gsk3=1) % Gsk3(axin=1),
     kf_axin_gsk3, kr_axin_gsk3)
Rule('axin_binds_apc', Axin(bcat=None, apc=None) + Apc(axin=None, aa15=None, aa20=None, state='u') |
     Axin(bcat=None, apc=1) % Apc(axin=1, aa15=None, aa20=None, state='u'),
     kf_axin_apc, kr_axin_apc)

# Bcat into dest comp
# Here, we are requiring that the aa20 site of APC is NOT bound to another beta-catenin molecule in order for
# beta-catenin to bind to Axin
# NOTE: we are only allowing Bcat to bind to Axin if APC is unphosphorylated (this completes the loop and allows
# the destruction complex to be recycled)
Rule('Bcat_binds_dtcpx',
     Bcat(top=None, bottom=None, nterm='u', state='x') +
     Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=1) % Apc(axin=1, aa20=None, aa15=None, state='u') |
     Bcat(top=2, bottom=3, nterm='u', state='x') %
     Axin(bcat=2, ck1a=ANY, gsk3=ANY, apc=1) % Apc(axin=1, aa20=None, aa15=3, state='u'),
     kf_bcat_dtcpx, kr_bcat_dtcpx)
# ###################

# STEP 2 ############
# NOTE: Bcat can only be bound to Axin at the 'top' site if ck1a, gsk3, and apc are all also bound
# Therefore, we don't need to explicitly include those three sites in the rule
# Also, dephosphorylation (reverse part of the rule) is implicitly modeled as due to PP2A
Rule('bcat_p_ck1a',
     Bcat(top=1, bottom=3, nterm='u') % Axin(bcat=1, apc=2) % Apc(axin=2, aa15=3, aa20=None, state='u') |
     Bcat(top=1, bottom=3, nterm='p1') % Axin(bcat=1, apc=2) % Apc(axin=2, aa15=3, aa20=None, state='u'),
     kf_bcat_phos_ck1a, k_dephos)  # ck1a phos first
# ###################

# STEP 3 #############
Rule('bcat_p_gsk3',
     Bcat(top=1, bottom=3, nterm='p1') % Axin(bcat=1, apc=2, gsk3=4) % Gsk3(axin=4, lithium=None) %
     Apc(axin=2, aa15=3, aa20=None, state='u') >>
     Bcat(top=1, bottom=3, nterm='p2') % Axin(bcat=1, apc=2, gsk3=4) % Gsk3(axin=4, lithium=None) %
     Apc(axin=2, aa15=3, aa20=None, state='u'),
     kf_bcat_phos_gsk3)  # then gsk3b phos

# Don't restrict dephosphorylation of Bcat to case where Gsk3 is not bound to lithium
# (i.e., Li could bind after phosphorylation of Bcat by Gsk3. Should be allowed to desphosphorylate in that case)
Rule('bcat_unp2',
     Bcat(top=1, bottom=3, nterm='p2') % Axin(bcat=1, apc=2) % Apc(axin=2, aa15=3, aa20=None, state='u') >>
     Bcat(top=1, bottom=3, nterm='p1') % Axin(bcat=1, apc=2) % Apc(axin=2, aa15=3, aa20=None, state='u'),
     k_dephos)
#######################

Rule('lithium_binds_GSK3', Gsk3(lithium=None) + Li(gsk3=None) | Gsk3(lithium=1) % Li(gsk3=1), kf_gsk3_li, kr_gsk3_li)

generate_equations(model)
for sp in model.species:
    print(sp)
# tspan = np.linspace(0, 5, 101)
# sim = ScipyOdeSimulator(model, tspan, verbose=False)
# x = sim.run()
# for i in range(len(model.species)):
#     plt.plot(tspan, x.all["__s%d" % i], lw=2, label="species %d" % i)
# plt.legend(loc=0)
# plt.show()
quit()

# TODO start from here
# Bcat binds to aa15 site of APC (it is assumed phosphorylation state of Bcat does not affect binding)
Rule('Bcat_binds_aa15', Bcat(top=1, bottom=None) % Axin(bcat=1, apc=2, gsk3=3, ck1a=4) %
     Apc(axin=2, state='u', aa15=None, aa20=None) % Gsk3(axin=3) % Ck1a(axin=4) | Bcat(top=1, bottom=5) %
     Axin(bcat=1, apc=2, gsk3=3, ck1a=4) % Apc(axin=2, state='u', aa15=5, aa20=None) % Gsk3(axin=3) % Ck1a(axin=4),
     kf_bcat_apc, kr_bcat_apc)


# APC phosphorylated by ck1a and gsk3
# Dephosphorylation is assumed to be due to PP2A (implicit)
Rule('apc_p_ck1a', Bcat(top=1, bottom=5, nterm='p2') % Axin(bcat=1, apc=2, gsk3=3, ck1a=4) %
     Gsk3(axin=3) % Ck1a(axin=4) % Apc(axin=2, aa15=5, state='u') >>
     Bcat(top=1, bottom=5, nterm='p2') % Axin(bcat=1, apc=2, gsk3=3, ck1a=4) %
     Gsk3(axin=3) % Ck1a(axin=4) % Apc(axin=2, aa15=5, state='p1'), kf_apc_phos_ck1a)


Rule('apc_unp1', Apc(state='p1', aa20=None) >> Apc(state='u', aa20=None), k_dephos)

Rule('apc_p_gsk3', Bcat(top=1, bottom=5, nterm='p2') % Apc(axin=2, aa15=5, state='p1') %
     Axin(bcat=1, apc=2, gsk3=3, ck1a=4) % Gsk3(axin=3, lithium=None) % Ck1a(axin=4) >>
     Bcat(top=1, bottom=5, nterm='p2') % Apc(axin=2, aa15=5, state='p2') % Axin(bcat=1, apc=2, gsk3=3, ck1a=4) %
     Gsk3(axin=3, lithium=None) % Ck1a(axin=4),
     kf_apc_phos_gsk3)

Rule('apc_unp2', Apc(state='p2', aa20=None) >> Apc(state='p1', aa20=None), k_dephos)

# phosphorylation forces bcat detach from axin

Rule('bcat_binds_apc', Bcat(top=1, nterm='p2', bottom=3) % Axin(bcat=1, apc=4) % Ck1a() % Gsk3() %
     Apc(axin=4, aa20=None, state='p2', aa15=3) >> Bcat(top=2, nterm='p2', bottom=None) % Axin(bcat=None, apc=4) % Ck1a() % Gsk3() %
     Apc(axin=4, aa20=2, state='p2', aa15=None), kf_bcat_binds_apc)




# apc is still bound to axin and bcat and apc are both phos

#  Btrcp binds bcat
# we are assuming that btrcp can bind (and unbind) to bcat whether it is ubiquitinated or not
Rule('btrcp_binds_bcat',
     Bcat(top=1, nterm='p2', bottom=None) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) + Btrcp(bcat=None) |
     Bcat(top=1, nterm='p2', bottom=2) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) % Btrcp(bcat=2),
     kf_btrcp_binds_bcat, kr_btrcp_binds_bcat)

# should axin=any be included?

#  Btrcp ubiquitinates bcat
Rule('Bcat_ubiq',
     Bcat(top=1, bottom=2, state='x') % Apc(aa20=1) % Btrcp(bcat=2) >>
     Bcat(top=1, bottom=2, state='ub') % Apc(aa20=1) % Btrcp(bcat=2), k_bcat_ubiq)

# Apc release bcat by dephos
Rule('apc_release_bcat',
     Bcat(top=1, bottom=2, state='ub') % Apc(aa20=1, state='p2') % Btrcp(bcat=2) >>
     Bcat(top=None, bottom=2, state='ub') % Btrcp(bcat=2) + Apc(aa20=None, state='u'), k_bcat_release)

# Bcat degraded by proteosome
Rule('bcat_degradation', Bcat(top=None, bottom=2, state='ub') % Btrcp(bcat=2) >> Btrcp(bcat=None), k_bcat_deg)

if __name__ == '__main__':

     #running simulations

     '''
     tspan = np.linspace(0, 1, 101)
     sim = ScipyOdeSimulator(model, tspan, verbose=True)
     result = sim.run()
     # for obs in model.observables:
          # if obs.name == 'bcat_free':
          #      plt.legend(loc=0)
          #      plt.figure()
     for i, group in enumerate(OBS):
         plt.figure()
         for obs_name in group:
             plt.plot(tspan, result.observables[obs_name], lw=2, label=obs_name)
         plt.legend(loc=0)
         plt.xlabel('time')
         plt.ylabel('concentration')
         ###
         if i == len(OBS)-1:
             plt.ylim(ymin=99.9875, ymax=100.0025)
     '''

     # In silico Li experiments of GSK3 activity
     # Comparing to Fig. 1(a) of Stambolic et al. (1996): doi:10.1016/s0960-9822(02)70790-2
     tspan = np.linspace(0, 40, 101)
     sim = ScipyOdeSimulator(model, tspan, verbose=False)
     # print('monomers %d' % len(model.monomers))
     # print('rules %d' % len(model.rules))
     # print('species %d' % len(model.species))
     # print('reactions %d' % lenw(model.reactions))
     names = []
     for rule in model.rules:
         names.append(rule.name)
     names = [n.lower() for n in names]
     names.sort()
     for name in names:
         print(name)
     quit()
     result = sim.run()

     # plt.plot(tspan, result.observables['GSK3_activity'], lw=2, label='GSK3_activity')
     #'kf_bcat_phos_gsk3'

     #'kf_gsk3_li' [0.1,1]

     Li_conc = np.arange(0, 101, 5)
     for kf in [1,10,100]:
         gsk3_activity = []
         for conc in Li_conc:
             print(conc, kf)
             result = sim.run(param_values={'Li_0': conc,'kf_bcat_phos_gsk3': kf})
             gsk3_activity.append(result.observables['GSK3_activity'][-1])

         gsk3_activity = np.array(gsk3_activity)

         plt.plot(Li_conc, gsk3_activity/gsk3_activity[0], 'o', label="kf=%g" % kf)
     plt.xlabel('Li concentration')
     plt.ylabel('GSK3 activity')
     plt.legend(loc=0)

     plt.show()


