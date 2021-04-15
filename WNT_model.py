from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator


Model()

# monomers

Monomer('Bcat', ['tcf4', 'gsk3b', 'apc', 'state', 'btrcp', 'loc'], {'state': ['x', 'ub'], 'loc' : ['cyt', 'nuc']})  # Bcat
Monomer('Gli2', ['btrcp', 'gPthlh', 'state', 'loc'], {'state': ['x', 'p', 'ub'], 'loc' : ['cyt', 'nuc']})  # Gli2 protein
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
Monomer('Apc', ['bcat'])   # apc
Monomer('Dkk1', ['rec'])   # wnt receptor inhibitor
Monomer('Wif1', ['wnt']) # wnt ligand inhibitor

# initials

Initial(Bcat(tcf4=None,apc=None,gsk3b=None,state='x',btrcp=None), Parameter('Bcat_0', 50))
Initial(Apc(bcat=None), Parameter('Apc_0', 500))
Initial(Gsk3b(bcat=None,dvl=None), Parameter('Gsk3b_0', 50))
Initial(Dvl(rec=1,gsk3b=None) % Rec(dvl=1,wnt=2) % Wnt(rec=2), Parameter('Dvl_rec_wnt_0', 50))
Initial(Gli2(btrcp=None,state='ub',gPthlh=None), Parameter('Gli2_0', 50))
Initial(Btrcp(b=None), Parameter('Btrcp_0', 50))

# Observables
Observable('pthrp_tot', Pthrp())
Observable('gli2_free', Gli2(btrcp=None,gPthlh=None))

# parameter

Parameter('kubiq', 15)
Parameter('knaked', 2)
Parameter('kfree', 47)
Parameter('kpromg',12)
Parameter('kgprod', 14)
Parameter('kINH', .22)
Parameter('k_beta_nuc', 1)

# rules

### ADD loc TO ALL B-CAT RULES ###

# Beta catenin ubiquitination and release from destruction complex
Rule('Bcat_u', Bcat(gsk3b=1,apc=2,state='x',btrcp=None) % Gsk3b(bcat=1) % Apc(bcat=2) + Btrcp(b=None) \
     >> Bcat(gsk3b=None,apc=None,state='ub',btrcp=1) % Btrcp(b=1) + Gsk3b(bcat=None) + Apc(bcat=None), kubiq)

# Beta catenin degradation
Rule('Bcat_degradation', Bcat(gsk3b=None,apc=None,state='ub',btrcp=1) % Btrcp(b=1) >> Btrcp(b=None), k_bcat_deg)

# Beta cantenin in destruction complex binds to DVL at the receptor
Rule('Bcat_naked', Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=None) % Apc(bcat=2) + Dvl(gsk3b=None,rec=ANY) >> \
     Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Dvl(gsk3b=3,rec=ANY), knaked)

# Release of APC from destruction complex by DVL
Rule('Bcat_apc', Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Dvl(gsk3b=3,rec=ANY) >> \
     Bcat(gsk3b=1,apc=None) % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) + Apc(bcat=None), kapc)

# Beta catenin released from destruction complex and DVL
Rule('Bcat_gsk3b', Bcat(gsk3b=1,apc=None) % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) >> Bcat(gsk3b=None,apc=None) \
     + Gsk3b(bcat=None,dvl=None) + Dvl(gsk3b=None,rec=ANY), kgsk)

# Beta catenin translocation to the nucleus
Rule('Bcat_to_nucleus', Bcat(gsk3b=None,apc=None,loc=cyt) >> Bcat(gsk3b=None,apc=None,loc=nuc), k_beta_nuc)

# Beta catenin in nuclues binds to TCF4
Rule('Bcat_tcf4', Bcat(tcf4=None,loc=nuc) + Tcf4(bcat=None) >> Bcat(tcf4=1,loc=nuc) % Tcf4(bcat=1), ktcf4)

# TCF4 binds to promoter of Gli2 gene
Rule('tcf4_binds_gGli2', Tcf4(gGli2=None) + gGli2(tcf4=None) | Tcf4(gGli2=1) % gGli2(tcf4=1), \
     kf_tcf4_gli2, kr_tcf4_gli2)

# SMAD3 binds to promoter of Gli2 gene
Rule('smad3_binds_gGli2', Smad3(gGli2=None) + gGli2(smad3=None) | Smad3(gGli2=1) % gGli2(smad3=1), \
     kf_smad3_gli2, kr_smad3_gli2)

# Gli2 gene transcription
Rule('gli2_transcription', gGli2(tcf4=1,smad3=ANY) % Tcf4(gGli2=1,bcat=ANY) >> \
     gGli2(tcf4=1,smad3=ANY) % Tcf4(gGli2=1,bcat=ANY) + mGli2(), k_gli2_tx)

# Gli2 translation
Rule('gli2_translation', mGli2() >> mGli2() + Gli2(btrcp=None, gPthlh=None, state='x', loc='cyt'), k_gli2_tl)

# TCF4 binds to promoter of PTHlH gene
Rule('tcf4_binds_gPthlh', Tcf4(gPthlh=None) + gPthlh(tcf4=None) | Tcf4(gPthlh=1) % gPthlh(tcf4=1), \
     kf_tcf4_pthlh, kr_tcf4_pthlh)

# Gli2 binds to promoter of PTHlH gene
Rule('tcf4_binds_gPthlh', Gli2(gPthlh=None,loc=nuc) + gPthlh(gli2=None) | Gli2(gPthlh=1,loc=nuc) % gPthlh(gli2=1), \
     kf_gli2_pthlh, kr_gli2_pthlh)

# PTHlH gene transcription
Rule('gPthlh', gPthlh(tcf4=1,gli2=ANY) % Tcf4(gPthlh=1,bcat=ANY) >> \
     >> gPthlh(tcf4=1,gli2=ANY) % Tcf4(gPthlh=1,bcat=ANY) + mPthlh(), k_pthlh_tx)

# PTHrP translation
Rule('Pthrp_translation', mPthlh() >> mPthlh() + Pthrp(), k_pthrp_tl)

# Gli2 phosporylation
Rule('Gli2_phospo', Gli2(state='x', loc='cyt', btrcp=None) >> Gli2(state='p', loc='cyt', btrcp=None), k_gli2_phos)

# Phopho-Gli2 translocates to nuclues
Rule('Gli2_to_nuc', Gli2(state='p', loc='cyt') >> Gli2(state='p', loc='nuc'), k_gli2_nuc)

# Gli2 binds BTRCP
Rule('Gli2_binds_BTRCP', Gli2(state='x', loc='cyt', btrcp=None) + Btrcp(gli2=None) | \
     Gli2(state='x', loc='cyt', btrcp=1) % Btrcp(gli2=1), kf_gli2_btrcp, kr_gli2_btrcp)

# Gli2 ubiquitination
Rule('Gli2_ubiq', Gli2(state='x', loc='cyt', btrcp=ANY) >> Gli2(state='ub', loc='cyt', btrcp=ANY), k_gli2_ubiq)

# Gli2 degradation
Rule('Gli2_degradation', Gli2(state='ub', loc='cyt', btrcp=1) % Btrcp(gli2=1) >> Btrcp(gli2=None), k_gli2_deg)

# WNT3A ligand binds receptor
Rule('wnt_binds_rec', Wnt(rec=None) + Rec(wnt=None) | Wnt(rec=1) % Rec(wnt=1), kf_wnt_rec, kr_wnt_rec)

# DVL binds WNT-bound receptor
Parameter('kf_dvl_rec', 10)
Parameter('kr_dvl_rec', 10)
Rule('dvl_binds_rec', Dvl(rec=None) + Rec(wnt=1,dvl=None) % Wnt(rec=1) | Dvl(rec=2) % Rec(wnt=1,dvl=2) % Wnt(rec=1) \
     kf_dvl_rec, kr_dvl_rec)

# DVL binds unbound receptor (only happens for high rigidity)
Parameter('kf_dvl_rec_rigid', 0)
Parameter('kr_dvl_rec_rigid', 0)
Rule('dvl_binds_rec', Dvl(rec=None) + Rec(wnt=None,dvl=None) | Dvl(rec=1) % Rec(wnt=None,dvl=1) \
     kf_dvl_rec_rigid, kr_dvl_rec_rigid)

# DKK1 binds receptor
Rule('dkk1_binds_rec', Dkk1(rec=None) + Rec(wnt=None) | Dkk1(rec=1) % Rec(wnt=1), kf_dkk1_rec, kr_dkk1_rec)

# WIF binds WNT3A
Rule('wif_binds_wnt', Wif1(wnt=None) + Wnt(rec=None) | Wif1(wnt=1) % Wnt(rec=1), kf_wif_wnt, kr_wif_wnt)



# from pysb.bng import generate_equations
# generate_equations(model,verbose=True)




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