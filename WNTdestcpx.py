from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator

Model()
# monomer
#destcpx monomers
Monomer('Axin', ['bcat', 'gsk3', 'ck1a', 'apc'])
Monomer('Gsk3', ['axin', 'lithium', 'dvl']) #added dvl component
Monomer('Ck1a', ['axin'])
Monomer('Apc', ['axin', 'aa15', 'aa20', 'state'], {'state': ['u', 'p1', 'p2']})
Monomer('Bcat', ['top', 'bottom', 'nterm', 'state', 'tcf4', 'loc'], #Added tcf4 and loc components
        {'state': ['x', 'ub'],'nterm': ['u', 'p1', 'p2'], 'loc' : ['cyt', 'nuc']}) # p1 by ck1a, p2 by gsk3
Monomer('Btrcp', ['bcat', 'b'])
Monomer('Li', ['gsk3'])
#WNT model monomers
Monomer('Gli2', ['btrcp', 'g_pthlh', 'state', 'loc'], {'state': ['x', 'p', 'ub'], 'loc' : ['cyt', 'nuc']})  # Gli2 protein
Monomer('gGli2', ['tcf4', 'smad3'])  # Gli2 gene
Monomer('mGli2')  # Gli2 mRNA
Monomer('Wnt', ['rec'])  # Wnt, rec is receptor complex - lrp5 and fzd
Monomer('Tcf4', ['g_gli2', 'g_pthlh', 'bcat'])  # Tcf4
Monomer('Smad3', ['g_gli2'])  # Smad3
Monomer('Rec', ['wnt', 'dvl'])  # Receptor
Monomer('Dvl', ['rec', 'gsk3'])  # Dvl
Monomer('gPthlh', ['gli2', 'tcf4'])  # Pthlh gene
Monomer('mPthlh')  # Pthlh mRNA
Monomer('Pthrp')  # Pthrp protein
Monomer('Dkk1', ['rec'])   # wnt receptor inhibitor
Monomer('Wif1', ['wnt']) # wnt ligand inhibitor

# Initials
Parameter('Bcat_0', 100)
Parameter('Gli2_0', 50)
Parameter('gGli2_0', 10)
Parameter('Wnt_0', 100)
Parameter('Btrcp_0', 50)
Parameter('Tcf4_0', 100)
Parameter('Smad3_0', 50)
Parameter('Gsk3_0', 50) #0
Parameter('Axin_0', 100) #0
Parameter('Rec_0', 100)
Parameter('Dvl_0', 50)
Parameter('gPthlh_0', 1)
Parameter('Apc_0', 50) #0
Parameter('Dkk1_0', 100)
Parameter('Wif1_0', 0)
Parameter('Ck1a_0', 50)
Parameter('Li_0', 0)

Initial(Axin(bcat=None, gsk3=None, ck1a=None, apc=None), Axin_0)
Initial(Gsk3(axin=None, lithium=None, dvl=None), Gsk3_0)
Initial(Ck1a(axin=None), Ck1a_0)
Initial(Apc(axin=None, aa15=None, aa20=None, state='u'), Apc_0)
Initial(Bcat(tcf4=None,top=None,bottom=None,state='x',loc='cyt',nterm='u'), Bcat_0)
Initial(Btrcp(bcat=None, b=None), Btrcp_0)
Initial(Li(gsk3=None), Li_0)
Initial(Gli2(btrcp=None,state='ub',g_pthlh=None,loc='cyt'), Gli2_0)
Initial(gGli2(tcf4=None,smad3=None), gGli2_0)
Initial(Wnt(rec=None), Wnt_0)
Initial(Tcf4(g_gli2=None, g_pthlh=None, bcat=None), Tcf4_0)
Initial(Smad3(g_gli2=None), Smad3_0)
Initial(Rec(wnt=None,dvl=None), Rec_0)
Initial(Dvl(rec=None,gsk3=None), Dvl_0)
Initial(gPthlh(gli2=None,tcf4=None), gPthlh_0)
#Initial(Dvl(rec=1,gsk3=None) % Rec(dvl=1,wnt=2) % Wnt(rec=2), Parameter('Dvl_rec_wnt_0', 50))
Initial(Dkk1(rec=None), Dkk1_0)
Initial(Wif1(wnt=None), Wif1_0)

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
Observable('apc_p1', Apc(state='p1'))
Observable('apc_p2', Apc(state='p2'))  # We think phos of Apc by GSK3 is "kinase activity" in Stambolic1996
Observable('bcat_apc_aa15', Apc(aa15=ANY, aa20=None))
Observable('bcat_apc_aa20', Apc(aa15=None, aa20=ANY))
Observable('bcat_apc_both', Apc(aa15=ANY, aa20=ANY))
Observable('bcat_ub_apc', Bcat(top=1,bottom=None,state='ub') % Apc(aa20=1))
Observable('bcat_ub_apc_btrcp', Bcat(top=1, bottom=2, state='ub') % Apc(aa20=1) % Btrcp(bcat=2))
Observable('bcat_ub_btrcp', Bcat(top=None, bottom=1, state='ub') % Btrcp(bcat=1))
Observable('bcat_ub_btrcp_total', Bcat(bottom=1, state='ub') % Btrcp(bcat=1))
Observable('bcat_ub_total', Bcat(state='ub'))
Observable('bcat_total', Bcat())
Observable('Li_total', Li())
Observable('Li_Gsk3', Li(gsk3=ANY))
###
Observable('GSK3_activity', Bcat(nterm='p2') + Apc(state='p2'))
destcpx_observables = [
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

Observable('pthrp_tot', Pthrp())
Observable('gli2_tot', Gli2())
Observable('gli2_cyt', Gli2(loc='cyt'))
Observable('gli2_nuc', Gli2(loc='nuc'))
Observable('bcat_tot', Bcat())
Observable('bcat_cyt', Bcat(loc='cyt'))
Observable('bcat_nuc', Bcat(loc='nuc'))
wnt_observables = [pthrp_tot,gli2_tot,gli2_cyt,gli2_nuc,bcat_tot,bcat_cyt,bcat_nuc]






# Rate constants
#destcpx parameters
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
Parameter('kf_apc_phos_gsk3', 200) #200
Parameter('k_dephos', 0.1)
Parameter('kf_bcat_binds_apc', 100)
Parameter('kf_btrcp_binds_bcat', 0.01)
Parameter('kr_btrcp_binds_bcat', 1)
Parameter('k_bcat_ubiq', 0.1) #1
Parameter('k_bcat_release', 10) #1
Parameter('k_bcat_deg', 0.001) #0.1
Parameter('kf_gsk3_li', 1)
Parameter('kr_gsk3_li', 0.01)
#WNT parameters
k_bcat_dvl = [
     Parameter('kf_bcat_dvl', 10),
     Parameter('kr_bcat_dvl', 10)]
Parameter('k_bcat_apc', 1)
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





# Rules
def destcpx_rules():
    # Axin binding rules
    # We require beta-catenin to NOT be bound for these binding event to occur
    # We assume beta-catenin can only bind when all three of Ck1a, Gsk3, and Apc are bound
    Rule('axin_binds_ck1a', Axin(bcat=None, ck1a=None) + Ck1a(axin=None) >> Axin(bcat=None, ck1a=1) % Ck1a(axin=1),
         kf_axin_ck1a)
    Rule('axin_binds_gsk3', Axin(bcat=None, gsk3=None) + Gsk3(axin=None, dvl=None) >> Axin(bcat=None, gsk3=1) % Gsk3(axin=1, dvl=None),
         kf_axin_gsk3)
    Rule('axin_binds_apc', Axin(bcat=None, apc=None) + Apc(axin=None) >> Axin(bcat=None, apc=1) % Apc(axin=1),
         kf_axin_apc)

    # Axin unbinding rules
    # We only allow Ck1a, Gsk3, and Apc to unbind if at least one of the others is not bound. If all three are
    # bound we assume the complex never breaks apart

    # Ck1a unbinds
    Rule('axin_unbinds_ck1a', Axin(bcat=None, ck1a=1, gsk3=None, apc=None) % Ck1a(axin=1) >>
         Axin(bcat=None, ck1a=None, gsk3=None, apc=None) + Ck1a(axin=None), kr_axin_ck1a)

    Rule('axin_gsk3_unbinds_ck1a', Axin(bcat=None, ck1a=1, gsk3=ANY, apc=None) % Ck1a(axin=1) >>
         Axin(bcat=None, ck1a=None, gsk3=ANY, apc=None) + Ck1a(axin=None), kr_axin_ck1a)

    Rule('axin_apc_unbinds_ck1a', Axin(bcat=None, ck1a=1, gsk3=None, apc=ANY) % Ck1a(axin=1) >>
         Axin(bcat=None, ck1a=None, gsk3=None, apc=ANY) + Ck1a(axin=None), kr_axin_ck1a)

    # Gsk3 unbinds
    Rule('axin_unbinds_gsk3', Axin(bcat=None, gsk3=1, ck1a=None, apc=None) % Gsk3(axin=1, dvl=None) >>
         Axin(bcat=None, gsk3=None, ck1a=None, apc=None) + Gsk3(axin=None,dvl=None), kr_axin_gsk3)

    Rule('axin_ck1a_unbinds_gsk3', Axin(bcat=None, gsk3=1, ck1a=ANY, apc=None) % Gsk3(axin=1,dvl=None) >>
         Axin(bcat=None, gsk3=None, ck1a=ANY, apc=None) + Gsk3(axin=None,dvl=None), kr_axin_gsk3)

    Rule('axin_apc_unbinds_gsk3', Axin(bcat=None, gsk3=1, ck1a=None, apc=ANY) % Gsk3(axin=1,dvl=None) >>
         Axin(bcat=None, gsk3=None, ck1a=None, apc=ANY) + Gsk3(axin=None,dvl=None), kr_axin_gsk3)

    # APC unbinds
    Rule('axin_unbinds_apc', Axin(bcat=None, apc=1, gsk3=None, ck1a=None) % Apc(axin=1) >>
         Axin(bcat=None, apc=None, gsk3=None, ck1a=None) + Apc(axin=None), kr_axin_apc)

    Rule('axin_gsk3_unbinds_apc', Axin(bcat=None, apc=1, gsk3=ANY, ck1a=None) % Apc(axin=1) >>
         Axin(bcat=None, apc=None, gsk3=ANY, ck1a=None) + Apc(axin=None), kr_axin_apc)

    Rule('axin_ck1a_unbinds_apc', Axin(bcat=None, apc=1, gsk3=None, ck1a=ANY) % Apc(axin=1) >>
         Axin(bcat=None, apc=None, gsk3=None, ck1a=ANY) + Apc(axin=None), kr_axin_apc)

    # Bcat into dest comp
    # Here, we are requiring that the aa20 site of APC is NOT bound to another beta-catenin molecule in order for
    # beta-catenin to bind to Axin
    # NOTE: we are only allowing Bcat to bind to Axin if APC is unphosphorylated (this completes the loop and allows
    # the destruction complex to be recycled)
    Rule('Bcat_binds_dtcpx',
         Bcat(top=None, bottom=None, tcf4=None, loc='cyt', nterm='u') +
         Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=1) % Apc(axin=1, aa20=None, state='u') |
         Bcat(top=2, bottom=None, tcf4=None, loc='cyt', nterm='u') % Axin(bcat=2, ck1a=ANY, gsk3=ANY, apc=1) % Apc(axin=1, aa20=None,
                                                                                             state='u'),
         kf_bcat_dtcpx, kr_bcat_dtcpx)

    # We think Bcat can be phosphorylated and dephosphorylated when bound to Axin, whether it's bound to APC or not

    # NOTE: Bcat can only be bound to Axin at the 'top' site if ck1a, gsk3, and apc are all also bound
    # Therefore, we don't need to explicitly include those three sites in the rule
    # Also, dephosphorylation (reverse part of the rule) is implicitly modeled as due to PP2A
    Rule('bcat_p_ck1a', Bcat(top=1, tcf4=None, loc='cyt', nterm='u') % Axin(bcat=1) | Bcat(top=1, tcf4=None, loc='cyt', nterm='p1') % Axin(bcat=1),
         kf_bcat_phos_ck1a, k_dephos)  # ck1a phos first

    Rule('bcat_p_gsk3', Bcat(top=1, tcf4=None, loc='cyt', nterm='p1') % Axin(bcat=1, gsk3=2) % Gsk3(axin=2, lithium=None, dvl=None) >>
         Bcat(top=1, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1, gsk3=2) % Gsk3(axin=2, lithium=None, dvl=None),
         kf_bcat_phos_gsk3)  # then gsk3b phos

    # Don't restrict dephosphorylation of Bcat to case where Gsk3 is not bound to lithium
    # (i.e., Li could bind after phosphorylation of Bcat by Gsk3. Should be allowed to desphosphorylate in that case)
    Rule('bcat_unp2', Bcat(top=1, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1, gsk3=2) % Gsk3(axin=2, dvl=None) >>
         Bcat(top=1, tcf4=None, loc='cyt', nterm='p1') % Axin(bcat=1, gsk3=2) % Gsk3(axin=2, dvl=None), k_dephos)

    Rule('lithium_binds_GSK3', Gsk3(lithium=None, dvl=None) + Li(gsk3=None) | Gsk3(lithium=1, dvl=None) % Li(gsk3=1), kf_gsk3_li,
         kr_gsk3_li)

    # Bcat binds to aa15 site of APC (it is assumed phosphorylation state of Bcat does not affect binding)
    Rule('Bcat_binds_aa15', Bcat(top=1, bottom=None, tcf4=None, loc='cyt') % Axin(bcat=1, apc=2) % Apc(axin=2, state='u', aa15=None) |
         Bcat(top=1, bottom=3, tcf4=None, loc='cyt') % Axin(bcat=1, apc=2) % Apc(axin=2, state='u', aa15=3), kf_bcat_apc, kr_bcat_apc)

    # APC phosphorylated by ck1a and gsk3
    # Dephosphorylation is assumed to be due to PP2A (implicit)
    Rule('apc_p_ck1a', Bcat(bottom=1, tcf4=None, loc='cyt', nterm='p2') % Apc(aa15=1, state='u', axin=ANY) % Ck1a() % Gsk3(dvl=None) >>
         Bcat(bottom=1, tcf4=None, loc='cyt', nterm='p2') % Apc(aa15=1, state='p1', axin=ANY) % Ck1a() % Gsk3(dvl=None), kf_apc_phos_ck1a)

    Rule('apc_unp1', Apc(state='p1', aa20=None) >> Apc(state='u', aa20=None), k_dephos)

    Rule('apc_p_gsk3', Bcat(bottom=1, tcf4=None, loc='cyt', nterm='p2') % Apc(aa15=1, state='p1', axin=ANY) % Ck1a() % Gsk3(lithium=None,dvl=None) >>
         Bcat(bottom=1, tcf4=None, loc='cyt', nterm='p2') % Apc(aa15=1, state='p2', axin=ANY) % Ck1a() % Gsk3(lithium=None, dvl=None), kf_apc_phos_gsk3)

    Rule('apc_unp2', Apc(state='p2', aa20=None) >> Apc(state='p1', aa20=None), k_dephos)

    # phosphorylation forces bcat detach from axin

    Rule('bcat_binds_apc', Bcat(top=1, tcf4=None, loc='cyt', nterm='p2', bottom=3) % Axin(bcat=1) % Ck1a() % Gsk3(dvl=None) %
         Apc(aa20=None, state='p2', aa15=3) >> Bcat(top=2, tcf4=None, loc='cyt', nterm='p2', bottom=None) % Axin(
        bcat=None) % Ck1a() % Gsk3(dvl=None) %
         Apc(aa20=2, state='p2', aa15=None), kf_bcat_binds_apc)

    # apc is still bound to axin and bcat and apc are both phos

    #  Btrcp binds bcat
    # we are assuming that btrcp can bind (and unbind) to bcat whether it is ubiquitinated or not
    Rule('btrcp_binds_bcat',
         Bcat(top=1, tcf4=None, loc='cyt', nterm='p2', bottom=None) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) + Btrcp(bcat=None) |
         Bcat(top=1, tcf4=None, loc='cyt', nterm='p2', bottom=2) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) % Btrcp(bcat=2),
         kf_btrcp_binds_bcat, kr_btrcp_binds_bcat)

    # should axin=any be included?

    #  Btrcp ubiquitinates bcat
    Rule('Bcat_ubiq',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='x') % Apc(aa20=1) % Btrcp(bcat=2) >>
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='ub') % Apc(aa20=1) % Btrcp(bcat=2), k_bcat_ubiq)

    # Apc release bcat by dephos
    Rule('apc_release_bcat',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='ub') % Apc(aa20=1) % Btrcp(bcat=2) >>
         Bcat(top=None, bottom=2, tcf4=None, loc='cyt', state='ub') % Btrcp(bcat=2) + Apc(aa20=None), k_bcat_release)

    # Bcat degraded by proteosome
    Rule('bcat_degradation', Bcat(top=None, bottom=2, tcf4=None, loc='cyt', state='ub') % Btrcp(bcat=2) >> Btrcp(bcat=None), k_bcat_deg)


def wntmodel_rules():

    Rule('Bcat_DVL', Bcat(top=1,bottom=2, tcf4=None, loc='cyt') % Gsk3(axin=3,dvl=None) % Apc(aa15=2) % Axin(bcat=1, gsk3=3) + Dvl(gsk3=None,rec=ANY) | \
         Bcat(top=1,bottom=2, tcf4=None, loc='cyt') % Gsk3(axin=3,dvl=4) % Apc(aa15=2) % Axin(bcat=1, gsk3=3) % Dvl(gsk3=4,rec=ANY), *k_bcat_dvl)

    # Release of APC from destruction complex by DVL
    Rule('Bcat_APC', Bcat(top=1,bottom=2, tcf4=None, loc='cyt') % Axin(bcat=1, gsk3=3) % Gsk3(axin=3, dvl=4) % Apc(aa15=2) % Dvl(gsk3=4,rec=ANY) >> \
         Bcat(top=1,bottom=None, tcf4=None, loc='cyt') % Axin(bcat=1, gsk3=3) % Gsk3(axin=3,dvl=4) % Dvl(gsk3=4,rec=ANY) + Apc(aa15=None), k_bcat_apc)

    # Beta catenin released from destruction complex and DVL
    Rule('Bcat_release', Bcat(top=1, tcf4=None, loc='cyt') % Axin(bcat=1, gsk3=3) % Gsk3(axin=3,dvl=4) % Dvl(gsk3=4,rec=ANY) >> \
         Bcat(top=None, tcf4=None, loc='cyt') + Axin(bcat=None, gsk3=None) + Gsk3(axin=None,dvl=None) + Dvl(gsk3=None,rec=ANY), k_bcat_release)

    # Beta catenin translocation to the nucleus
    Rule('Bcat_to_nucleus', Bcat(top=None,bottom=None,loc='cyt') >> Bcat(top=None,bottom=None,loc='nuc'), \
         k_bcat_nuc)


    # Beta catenin in nuclues binds to TCF4
    # TODO: Does tcf4 bind gGli2 before bcat binds?
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
    Rule('Gli2_ubiq', Gli2(state='x', loc='cyt', btrcp=ANY) >> Gli2(state='ub', loc='cyt', btrcp=ANY), k_gli2_ubiq)

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

destcpx_rules()
wntmodel_rules()

#run simulation
tspan=np.linspace(0,40,101)
sim=ScipyOdeSimulator(model,tspan,verbose=False)
result=sim.run()

'''
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
'''

# run simulation(WNT)
# tspan=np.linspace(0,500,501)
# sim=ScipyOdeSimulator(model,tspan,verbose=True)
# traj=sim.run()
fig, axs=plt.subplots(nrows=4,ncols=2,figsize=(6.4,9.6))
row=0
col=0
for obs in wnt_observables:
     #plt.figure()
     axs[row,col].plot(tspan, result.observables[obs.name], lw=2, label=obs.name)
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






