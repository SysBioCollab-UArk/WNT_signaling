from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator

Model()

# monomers
Monomer('Bcat', ['tcf4', 'gsk3b', 'apc', 'state', 'btrcp', 'loc','axin'], {'state': ['x', 'ub'], 'loc' : ['cyt', 'nuc']})  # Bcat
Monomer('Gli2', ['btrcp', 'g_pthlh', 'state', 'loc'], {'state': ['x', 'p', 'ub'], 'loc' : ['cyt', 'nuc']})  # Gli2 protein
Monomer('gGli2', ['tcf4', 'smad3'])  # Gli2 gene
Monomer('mGli2')  # Gli2 mRNA
Monomer('Wnt', ['rec'])  # Wnt, rec is receptor complex - lrp5 and fzd
Monomer('Btrcp', ['b'])  # Btrcp it binds to the Ub Bcat
Monomer('Tcf4', ['g_gli2', 'g_pthlh', 'bcat'])  # Tcf4
Monomer('Smad3', ['g_gli2'])  # Smad3
Monomer('Gsk3b', ['bcat', 'dvl'])  # represents gsk3, axin, ck1
#we want to add rules for AXIN and CK1 so me also need monomers
Monomer('Axin',['bcat'])#
Monomer('Rec', ['wnt', 'dvl'])  # Receptor
Monomer('Dvl', ['rec', 'gsk3b'])  # Dvl
Monomer('gPthlh', ['gli2', 'tcf4'])  # Pthlh gene
Monomer('mPthlh')  # Pthlh mRNA
Monomer('Pthrp')  # Pthrp protein
Monomer('Apc', ['bcat'])   # apc
Monomer('Dkk1', ['rec'])   # wnt receptor inhibitor
Monomer('Wif1', ['wnt']) # wnt ligand inhibitor
#Monomer('LiCl',['Gsk3b'])
# initials
Parameter('Bcat_0', 100)
Parameter('Gli2_0', 50)
Parameter('gGli2_0', 1)
Parameter('Wnt_0', 100)
Parameter('Btrcp_0', 30) #50
Parameter('Tcf4_0', 100)
Parameter('Smad3_0', 50)
Parameter('Gsk3b_0', 50)
Parameter('Axin_0',100)
Parameter('Rec_0', 100)
Parameter('Dvl_0', 50)
Parameter('gPthlh_0', 1)
Parameter('Apc_0', 50)
Parameter('Dkk1_0', 100)
Parameter('Wif1_0', 0)

Initial(Bcat(apc=1,gsk3b=2,btrcp=None,axin=4,tcf4=None,state='x',loc='cyt') % Apc(bcat=1) % Gsk3b(bcat=2,dvl=None) % Axin(bcat=4), Bcat_0)
Initial(Gli2(btrcp=None,state='ub',g_pthlh=None,loc='cyt'), Gli2_0)
Initial(gGli2(tcf4=None,smad3=None), gGli2_0)
Initial(Wnt(rec=None), Wnt_0)
Initial(Btrcp(b=None), Btrcp_0)
Initial(Tcf4(g_gli2=None, g_pthlh=None, bcat=None), Tcf4_0)
Initial(Smad3(g_gli2=None), Smad3_0)
Initial(Gsk3b(bcat=None,dvl=None), Gsk3b_0)
Initial(Axin(bcat=None),Axin_0)
Initial(Rec(wnt=None,dvl=None), Rec_0)
Initial(Dvl(rec=None,gsk3b=None), Dvl_0)
Initial(gPthlh(gli2=None,tcf4=None), gPthlh_0)
Initial(Apc(bcat=None), Apc_0)
# Initial(Dvl(rec=1,gsk3b=None) % Rec(dvl=1,wnt=2) % Wnt(rec=2), Parameter('Dvl_rec_wnt_0', 50))
Initial(Dkk1(rec=None), Dkk1_0)
Initial(Wif1(wnt=None), Wif1_0)

# observables
Observable('pthrp_tot', Pthrp())
Observable('gli2_tot', Gli2())
Observable('gli2_cyt', Gli2(loc='cyt'))
Observable('gli2_nuc', Gli2(loc='nuc'))
Observable('bcat_tot', Bcat())
Observable('bcat_cyt', Bcat(loc='cyt'))
Observable('bcat_nuc', Bcat(loc='nuc'))

# rate constants
Parameter('k_bcat_ubiq', 0.1)
Parameter('k_bcat_deg', 0.001)
k_bcat_dvl = [
     Parameter('kf_bcat_dvl', 10),
     Parameter('kr_bcat_dvl', 10)]
Parameter('k_bcat_apc', 1)
Parameter('k_bcat_release', 10)
Parameter('k_bcat_nuc', 10)
k_bcat_tcf4 = [
     Parameter('kf_bcat_tcf4', 1),
     Parameter('kr_bcat_tcf4', 1)]
k_tcf4_gli2 = [
     Parameter('kf_tcf4_gli2', 10),
     Parameter('kr_tcf4_gli2', 1)]
k_smad3_gli2 = [
     Parameter('kf_smad3_gli2', 10),
     Parameter('kr_smad3_gli2', 1)]
Parameter('k_gli2_tx', 10)
Parameter('k_gli2_tl', 10)
Parameter('k_mgli2_deg', 1)
k_tcf4_pthlh = [
     Parameter('kf_tcf4_pthlh', 10),
     Parameter('kr_tcf4_pthlh', 1)]
k_gli2_pthlh = [
     Parameter('kf_gli2_pthlh', 10),
     Parameter('kr_gli2_pthlh', 1)]
Parameter('k_pthlh_tx', 10)
Parameter('k_pthrp_tl', 10)
Parameter('k_mPthlh_deg', 0.1)
Parameter('k_pthrp_deg', 0.1)
k_gli2_phos = [Parameter('kf_gli2_phos', 10),
               Parameter('kr_gli2_phos', 1)]
k_gli2_nuc = [Parameter('kf_gli2_nuc', 100),
              Parameter('kr_gli2_nuc', 1)]
k_gli2_btrcp = [
     Parameter('kf_gli2_btrcp', 10),
     Parameter('kr_gli2_btrcp', 1)]
Parameter('k_gli2_ubiq', 10)
Parameter('k_gli2_deg', 100)
k_wnt_rec = [
     Parameter('kf_wnt_rec', 10),
     Parameter('kr_wnt_rec', 1)]
k_dvl_rec = [
     Parameter('kf_dvl_rec', 10),
     Parameter('kr_dvl_rec', 1)]
k_dvl_rec_rigid = [
     Parameter('kf_dvl_rec_rigid', 100),
     Parameter('kr_dvl_rec_rigid', 1)]
k_dkk1_rec = [
     Parameter('kf_dkk1_rec', 1),
     Parameter('kr_dkk1_rec', 1)]
k_wif_wnt = [
     Parameter('kf_wif_wnt', 1),
     Parameter('kr_wif_wnt', 1)]

# rules

# Beta catenin ubiquitination and release from destruction complex (in the cytoplasm)
Rule('Bcat_ubiq', Bcat(gsk3b=1,apc=2,btrcp=None,state='x',loc='cyt') % Gsk3b(bcat=1) % Apc(bcat=2) + Btrcp(b=None) >> \
     Bcat(gsk3b=None,apc=None,btrcp=3,state='ub',loc='cyt') % Btrcp(b=3) + Gsk3b(bcat=None) + Apc(bcat=None), \
     k_bcat_ubiq)

# Beta catenin degradation
Rule('Bcat_degradation', Bcat(gsk3b=None,apc=None,btrcp=1,state='ub') % Btrcp(b=1) >> Btrcp(b=None), k_bcat_deg)

# Beta cantenin in destruction complex binds to DVL at the receptor
#Adding gsk3b and AXIN and find out what happen with it

Rule('Bcat_DVL', Bcat(gsk3b=1,apc=2,axin=4) % Gsk3b(bcat=1,dvl=None) % Apc(bcat=2) % Axin(bcat=4) + Dvl(gsk3b=None,rec=ANY) | \
     Bcat(gsk3b=1,apc=2, axin=4) % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Axin(bcat=4) % Dvl(gsk3b=3,rec=ANY), *k_bcat_dvl)

# Release of APC from destruction complex by DVL
Rule('Bcat_APC', Bcat(gsk3b=1,apc=2) % Gsk3b(bcat=1,dvl=3) % Apc(bcat=2) % Dvl(gsk3b=3,rec=ANY) >> \
     Bcat(gsk3b=1,apc=None) % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) + Apc(bcat=None), k_bcat_apc)

# Beta catenin released from destruction complex and DVL
Rule('Bcat_release', Bcat(gsk3b=1,apc=None) % Gsk3b(bcat=1,dvl=3) % Dvl(gsk3b=3,rec=ANY) >> \
     Bcat(gsk3b=None,apc=None) + Gsk3b(bcat=None,dvl=None) + Dvl(gsk3b=None,rec=ANY), k_bcat_release)

# Beta catenin translocation to the nucleus
Rule('Bcat_to_nucleus', Bcat(gsk3b=None,apc=None,btrcp=None,loc='cyt') >> Bcat(gsk3b=None,apc=None,btrcp=None,loc='nuc'), \
     k_bcat_nuc)

# Beta catenin in nuclues binds to TCF4
Rule('Bcat_tcf4', Bcat(tcf4=None,loc='nuc') + Tcf4(bcat=None) | Bcat(tcf4=1,loc='nuc') % Tcf4(bcat=1), *k_bcat_tcf4)

# TCF4 binds to promoter of Gli2 gene
Rule('tcf4_binds_gGli2', Tcf4(g_gli2=None) + gGli2(tcf4=None) | Tcf4(g_gli2=1) % gGli2(tcf4=1), *k_tcf4_gli2)

# SMAD3 binds to promoter of Gli2 gene
Rule('smad3_binds_gGli2', Smad3(g_gli2=None) + gGli2(smad3=None) | Smad3(g_gli2=1) % gGli2(smad3=1), *k_smad3_gli2)

# Gli2 gene transcription
Rule('gli2_transcription', gGli2(tcf4=1,smad3=ANY) % Tcf4(g_gli2=1,bcat=ANY) >> \
     gGli2(tcf4=1,smad3=ANY) % Tcf4(g_gli2=1,bcat=ANY) + mGli2(), k_gli2_tx)

# Gli2 translation
Rule('gli2_translation', mGli2() >> mGli2() + Gli2(btrcp=None, g_pthlh=None, state='x', loc='cyt'), k_gli2_tl)

# mGli2 degradation
Rule('mGli2_deg', mGli2() >> None, k_mgli2_deg)

# TCF4 binds to promoter of PTHlH gene
Rule('tcf4_binds_gPthlh', Tcf4(g_pthlh=None) + gPthlh(tcf4=None) | Tcf4(g_pthlh=1) % gPthlh(tcf4=1), *k_tcf4_pthlh)

# Gli2 binds to promoter of PTHlH gene
Rule('gli2_binds_gPthlh', Gli2(g_pthlh=None,loc='nuc') + gPthlh(gli2=None) | \
     Gli2(g_pthlh=1,loc='nuc') % gPthlh(gli2=1), *k_gli2_pthlh)

# PTHlH gene transcription
Rule('gPthlh_transcription', gPthlh(tcf4=1,gli2=ANY) % Tcf4(g_pthlh=1,bcat=ANY) >> \
     gPthlh(tcf4=1,gli2=ANY) % Tcf4(g_pthlh=1,bcat=ANY) + mPthlh(), k_pthlh_tx)

# PTHrP translation
Rule('Pthrp_translation', mPthlh() >> mPthlh() + Pthrp(), k_pthrp_tl)

# mPthlh degradation
Rule('mPthlh_deg', mPthlh() >> None, k_mPthlh_deg)

# PTHrP degrdation
Rule('Pthrp_deg', Pthrp() >> None, k_pthrp_deg)

# Gli2 phosporylation
Rule('Gli2_phospo', Gli2(state='x', loc='cyt', btrcp=None) | Gli2(state='p', loc='cyt', btrcp=None), *k_gli2_phos)

# Phopho-Gli2 translocates to nucleus
Rule('Gli2_to_nuc', Gli2(state='p', loc='cyt') | Gli2(state='p', loc='nuc'), *k_gli2_nuc)

# Gli2 binds BTRCP
Rule('Gli2_binds_BTRCP', Gli2(state='x', loc='cyt', btrcp=None) + Btrcp(b=None) | \
     Gli2(state='x', loc='cyt', btrcp=1) % Btrcp(b=1), *k_gli2_btrcp)

# Gli2 ubiquitination
Rule('Gli2_ubiq', Gli2(state='x', btrcp=ANY) >> Gli2(state='ub', btrcp=ANY), k_gli2_ubiq)

# Gli2 degradation
Rule('Gli2_degradation', Gli2(state='ub', btrcp=1) % Btrcp(b=1) >> Btrcp(b=None), k_gli2_deg)

# WNT3A ligand binds receptor
Rule('wnt_binds_rec', Wnt(rec=None) + Rec(wnt=None) | Wnt(rec=1) % Rec(wnt=1), *k_wnt_rec)

# DVL binds WNT-bound receptor
Rule('dvl_binds_wnt_rec', Dvl(rec=None) + Rec(wnt=1,dvl=None) % Wnt(rec=1) | \
     Dvl(rec=2) % Rec(wnt=1,dvl=2) % Wnt(rec=1), *k_dvl_rec)

# DVL binds unbound receptor (only happens for high rigidity)
Rule('dvl_binds_rec', Dvl(rec=None) + Rec(wnt=None,dvl=None) | Dvl(rec=1) % Rec(wnt=None,dvl=1), *k_dvl_rec_rigid)

# DKK1 binds receptor
Rule('dkk1_binds_rec', Dkk1(rec=None) + Rec(wnt=None) | Dkk1(rec=1) % Rec(wnt=1), *k_dkk1_rec)

# WIF binds WNT3A
Rule('wif_binds_wnt', Wif1(wnt=None) + Wnt(rec=None) | Wif1(wnt=1) % Wnt(rec=1), *k_wif_wnt)

# run simulation
tspan=np.linspace(0,40,101)
sim=ScipyOdeSimulator(model,tspan,verbose=True)
traj=sim.run()
fig, axs=plt.subplots(nrows=4,ncols=2,figsize=(6.4,9.6))
row=0
col=0
for obs in model.observables:
     #plt.figure()
     axs[row,col].plot(tspan, traj.observables[obs.name], lw=2, label=obs.name)
     axs[row,col].legend(loc=0)
     axs[row,col].set_xlabel('Time (arbitrary units)')
     axs[row,col].set_ylabel('Molecule count')
     axs[row,col].ticklabel_format(style='scientific')
     if ((col+1) % 2 == 0):
         row += 1
         col = 0
     else:
         col += 1

plt.tight_layout(pad=1)
plt.show()