from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator

Model()

# destcpx monomers
Monomer('Axin', ['bcat', 'gsk3', 'ck1a', 'apc'])
Monomer('Gsk3', ['axin', 'lithium', 'dvl'])  # added dvl component
Monomer('Ck1a', ['axin'])
Monomer('Apc', ['axin', 'aa15', 'aa20', 'state'], {'state': ['u', 'p1', 'p2']})
Monomer('Bcat', ['top', 'bottom', 'nterm', 'state', 'tcf4', 'loc'],  # added tcf4 and loc components
        {'state': ['x', 'ub'], 'nterm': ['u', 'p1', 'p2'], 'loc': ['cyt', 'nuc']})  # p1 by ck1a, p2 by gsk3
Monomer('Btrcp', ['b'])
Monomer('Li', ['gsk3'])

# WNT model monomers
Monomer('Rec', ['wnt', 'dvl'])
Monomer('Dvl', ['rec', 'gsk3'])
Monomer('Wnt', ['rec'])  # Wnt, rec is receptor complex - lrp5 and fzd
Monomer('Dkk1', ['rec'])  # Wnt receptor inhibitor
Monomer('Wif1', ['wnt'])  # Wnt ligand inhibitor
Monomer('Gli2', ['btrcp', 'g_pthlh', 'state', 'loc'], {'state': ['x', 'p', 'ub'], 'loc': ['cyt', 'nuc']})
Monomer('gGli2', ['tcf4', 'smad3'])  # Gli2 gene
Monomer('mGli2')  # Gli2 mRNA
Monomer('gPthlh', ['gli2', 'tcf4'])  # Pthlh gene
Monomer('mPthlh')  # Pthlh mRNA
Monomer('Pthrp')  # Pthrp protein
Monomer('Tcf4', ['g', 'bcat'])
Monomer('Smad3', ['g_gli2'])

# Initials
Parameter('Axin_0', 100)
Parameter('Gsk3_0', 50)
Parameter('Ck1a_0', 50)
Parameter('Apc_0', 50)
Parameter('Bcat_0', 0)  # 100)
Parameter('Btrcp_0', 1000)  # 30)
Parameter('Li_0', 0)

Parameter('Rec_0', 100)
Parameter('Dvl_0', 50)
Parameter('Wnt_0', 50)  # 50 todo
# Parameter('Dkk1_0', 100)
# Parameter('Wif1_0', 0)
Parameter('Gli2_0', 0)  # 50)
Parameter('gGli2_0', 1)
Parameter('gPthlh_0', 1)
Parameter('Tcf4_0', 100)
Parameter('Smad3_0', 50)

Initial(Axin(bcat=None, gsk3=None, ck1a=None, apc=None), Axin_0)
Initial(Gsk3(axin=None, lithium=None, dvl=None), Gsk3_0)
Initial(Ck1a(axin=None), Ck1a_0)
Initial(Apc(axin=None, aa15=None, aa20=None, state='u'), Apc_0)
Initial(Bcat(tcf4=None, top=None, bottom=None, state='x', loc='cyt', nterm='u'), Bcat_0)
Initial(Btrcp(b=None), Btrcp_0)
Initial(Li(gsk3=None), Li_0)

Initial(Rec(wnt=None, dvl=None), Rec_0)
Initial(Dvl(rec=None, gsk3=None), Dvl_0)
Initial(Wnt(rec=None), Wnt_0)
# Initial(Dkk1(rec=None), Dkk1_0)
# Initial(Wif1(wnt=None), Wif1_0)
Initial(Gli2(btrcp=None, state='x', g_pthlh=None, loc='cyt'), Gli2_0)
Initial(gGli2(tcf4=None, smad3=None), gGli2_0)
Initial(gPthlh(gli2=None, tcf4=None), gPthlh_0)
Initial(Tcf4(g=None, bcat=None), Tcf4_0)
Initial(Smad3(g_gli2=None), Smad3_0)

# Observables
Observable('Axin_free', Axin(bcat=None, gsk3=None, ck1a=None, apc=None))
Observable('Gsk3_free', Gsk3(axin=None, lithium=None, dvl=None))
Observable('Ck1a_free', Ck1a(axin=None))
Observable('Apc_free', Apc(axin=None, aa15=None, aa20=None, state='u'))
Observable('DestCpx', Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=ANY))



# Observable('Axin_tot', Axin())
# Observable('Axin_Ck1a', Axin(ck1a=ANY))
# Observable('Axin_Gsk3', Axin(gsk3=ANY))
# Observable('Axin_Apc', Axin(apc=ANY))
# Observable('Gsk3_tot', Gsk3())
# Observable('Ck1a_tot', Ck1a())
# Observable('Apc_u', Apc(state='u'))
# Observable('Apc_p1', Apc(state='p1'))
Observable('Apc_p2', Apc(state='p2'))  # We think phos of Apc by GSK3 is "kinase activity" in Stambolic 1996
Observable('Bcat_cyt', Bcat(loc='cyt'))
# Observable('Bcat_cyt_free', Bcat(top=None, bottom=None, loc='cyt'))
Observable('Bcat_nuc', Bcat(loc='nuc'))
# Observable('Bcat_Axin', Bcat(top=1) % Axin(bcat=1))
# Observable('Bcat_Axin_Apc', Bcat(top=1, bottom=2) % Axin(bcat=1) % Apc(aa15=2))
# Observable('Bcat_u', Bcat(nterm='u'))
# Observable('Bcat_p1', Bcat(nterm='p1'))
Observable('Bcat_p2', Bcat(nterm='p2'))
# Observable('Bcat_Apc_aa15', Apc(aa15=ANY, aa20=None))
# Observable('Bcat_Apc_aa20', Bcat(top=1) % Apc(aa20=1))
# Observable('Bcat_Apc_both', Apc(aa15=ANY, aa20=ANY))  # this should be zero
Observable('Bcat_x_Apc_aa20', Bcat(top=1, state='x') % Apc(aa20=1))
Observable('Bcat_ub_Apc_aa20', Bcat(top=1, state='ub') % Apc(aa20=1))
Observable('Bcat_ub_free', Bcat(top=None, bottom=ANY, state='ub'))
# Observable('Bcat_ub', Bcat(state='ub'))
Observable('Bcat_tot', Bcat())
Observable('DestCpx_Bcat', Axin(bcat=ANY, ck1a=ANY, gsk3=ANY, apc=ANY))
Observable('GSK3_activity', Bcat(nterm='p2') + Apc(state='p2'))
# Observable('Btrcp_tot', Btrcp())
# Observable('Li_total', Li())
# Observable('Li_Gsk3', Li(gsk3=ANY))

Observable('Gli2_tot', Gli2())
Observable('Gli2_cyt', Gli2(loc='cyt'))
Observable('Gli2_nuc', Gli2(loc='nuc'))
Observable('Pthrp_tot', Pthrp())

obs_to_plot = [
    ['Axin_free', 'Gsk3_free', 'Ck1a_free', 'Apc_free'],
    # ['Pthrp_tot']#,
    # ['Gli2_tot', 'Gli2_cyt', 'Gli2_nuc'],
    # ['Bcat_tot', 'Bcat_cyt', 'Bcat_nuc'],  # 'Bcat_cyt_free'],
    # ['Bcat_x_Apc_aa20', 'Bcat_ub_Apc_aa20', 'Bcat_ub_free'],
    # ['Apc_p2', 'Bcat_p2', 'GSK3_activity'],
    ['DestCpx', 'DestCpx_Bcat']
]

# Rate constants

# common parameters
Parameter('kf_btrcp_binds_bcat', 0.01)
Parameter('kr_btrcp_binds_bcat', 1)
Parameter('k_bcat_ubiq', 1)  # 1
Parameter('k_bcat_deg', 10)  # 0.1
Parameter('k_bcat_synth', 0.1)  # since beta-catenin is degraded by the proteosome, we need to add a source term

# destcpx parameters
Parameter('kf_axin_ck1a', 1)
Parameter('kr_axin_ck1a', 1e4)
Parameter('kf_axin_gsk3', 1)
Parameter('kr_axin_gsk3', 1e4)
Parameter('kf_axin_apc', 1)
Parameter('kr_axin_apc', 1e4)
Parameter('kf_bcat_dtcpx', 1)  # 100)
Parameter('kr_bcat_dtcpx', 100)  # 0.1)
# Parameter('kf_bcat_apc', 100)  # This parameter isn't used anymore. Bcat binds Axin and Apc at the same time now.
# Parameter('kr_bcat_apc', 0.1)  # This parameter isn't used anymore. Bcat binds Axin and Apc at the same time now.
Parameter('k_bcat_phos_ck1a', 10)
Parameter('k_bcat_phos_gsk3', 100)
Parameter('k_apc_phos_ck1a', 20)
Parameter('k_apc_phos_gsk3', 200)  # 200
Parameter('k_dephos', 0.1)
Parameter('k_bcat_binds_apc_aa20', 0.01)  # 0.005) todo
Parameter('k_bcat_ub_release_apc', 1)  # 10
Parameter('kf_gsk3_li', 1)
Parameter('kr_gsk3_li', 0.01)

# WNT parameters
k_bcat_dvl = [
     Parameter('kf_bcat_dvl', 1),  # 1), todo
     Parameter('kr_bcat_dvl', 1)]
Parameter('k_apc_release_dvl', 1)
Parameter('k_bcat_release_dvl', 1)
k_bcat_nuc = [
    Parameter('kf_bcat_nuc', 1),
    Parameter('kr_bcat_nuc', 1)]  # 1 todo
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
k_gli2_nuc = [Parameter('kf_gli2_nuc', 1),  # 100
              Parameter('kr_gli2_nuc', 10)]
k_gli2_btrcp = [
     Parameter('kf_gli2_btrcp', 10),
     Parameter('kr_gli2_btrcp', 1)]
Parameter('k_gli2_ubiq', 10)
Parameter('k_gli2_deg', 100)
k_wnt_rec = [
     Parameter('kf_wnt_rec', 1),
     Parameter('kr_wnt_rec', 10)]
k_dvl_rec = [
     Parameter('kf_dvl_rec', 1),
     Parameter('kr_dvl_rec', 10)]
k_dvl_rec_rigid = [
     Parameter('kf_dvl_rec_rigid', 0),  # 100),
     Parameter('kr_dvl_rec_rigid', 0)]  # 1)]
k_dkk1_rec = [
     Parameter('kf_dkk1_rec', 1),
     Parameter('kr_dkk1_rec', 1)]
k_wif_wnt = [
     Parameter('kf_wif_wnt', 1),
     Parameter('kr_wif_wnt', 1)]

# Common Rules

# Btrcp binds to Bcat after Bcat has detached from Axin
Rule('Btrcp_binds_Bcat',
     Bcat(top=1, bottom=None, tcf4=None, loc='cyt', nterm='p2', state='x') %
     Apc(aa15=None, aa20=1, state='p2', axin=ANY) + Btrcp(b=None) |
     Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p2', state='x') %
     Apc(aa15=None, aa20=1, state='p2', axin=ANY) % Btrcp(b=2),
     kf_btrcp_binds_bcat, kr_btrcp_binds_bcat)

# Beta catenin ubiquitination
Rule('Bcat_ubiqitination',
     Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='x') % Apc(aa20=1) % Btrcp(b=2) >>
     Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='ub') % Apc(aa20=1) % Btrcp(b=2),
     k_bcat_ubiq)

# Bcat degraded by proteosome
Rule('Bcat_degradation',
     Bcat(top=None, bottom=2, tcf4=None, loc='cyt', state='ub') % Btrcp(b=2) >> Btrcp(b=None),
     k_bcat_deg)

# Bcat synthesis
Rule('Bcat_synthesis', None >> Bcat(top=None, bottom=None, nterm='u', state='x', tcf4=None, loc='cyt'), k_bcat_synth)


def create_destcpx_rules(create=True):
    if not create:
        return False

    # Axin binding rules
    # We require beta-catenin to NOT be bound for these binding event to occur
    # We assume beta-catenin can only bind when all three of Ck1a, Gsk3, and Apc are bound
    Rule('Axin_binds_Ck1a',
         Axin(bcat=None, ck1a=None) + Ck1a(axin=None) >>
         Axin(bcat=None, ck1a=1) % Ck1a(axin=1),
         kf_axin_ck1a)

    Rule('Axin_binds_Gsk3',
         Axin(bcat=None, gsk3=None) + Gsk3(axin=None, dvl=None) >>
         Axin(bcat=None, gsk3=1) % Gsk3(axin=1, dvl=None),
         kf_axin_gsk3)

    Rule('Axin_binds_Apc',
         Axin(bcat=None, apc=None) + Apc(axin=None) >>
         Axin(bcat=None, apc=1) % Apc(axin=1),
         kf_axin_apc)

    # Axin unbinding rules
    # We only allow Ck1a, Gsk3, and Apc to unbind if at least one of the others is not bound. If all three are
    # bound we assume the complex never breaks apart

    # Ck1a unbinds
    Rule('Axin_unbinds_Ck1a',
         Axin(bcat=None, ck1a=1, gsk3=None, apc=None) % Ck1a(axin=1) >>
         Axin(bcat=None, ck1a=None, gsk3=None, apc=None) + Ck1a(axin=None),
         kr_axin_ck1a)

    Rule('Axin_Gsk3_unbinds_Ck1a',
         Axin(bcat=None, ck1a=1, gsk3=ANY, apc=None) % Ck1a(axin=1) >>
         Axin(bcat=None, ck1a=None, gsk3=ANY, apc=None) + Ck1a(axin=None),
         kr_axin_ck1a)

    Rule('Axin_Apc_unbinds_Ck1a',
         Axin(bcat=None, ck1a=1, gsk3=None, apc=ANY) % Ck1a(axin=1) >>
         Axin(bcat=None, ck1a=None, gsk3=None, apc=ANY) + Ck1a(axin=None),
         kr_axin_ck1a)

    # Gsk3 unbinds
    Rule('Axin_unbinds_Gsk3',
         Axin(bcat=None, gsk3=1, ck1a=None, apc=None) % Gsk3(axin=1, dvl=None) >>
         Axin(bcat=None, gsk3=None, ck1a=None, apc=None) + Gsk3(axin=None,dvl=None),
         kr_axin_gsk3)

    Rule('Axin_Ck1a_unbinds_Gsk3',
         Axin(bcat=None, gsk3=1, ck1a=ANY, apc=None) % Gsk3(axin=1,dvl=None) >>
         Axin(bcat=None, gsk3=None, ck1a=ANY, apc=None) + Gsk3(axin=None,dvl=None),
         kr_axin_gsk3)

    Rule('Axin_Apc_unbinds_Gsk3',
         Axin(bcat=None, gsk3=1, ck1a=None, apc=ANY) % Gsk3(axin=1,dvl=None) >>
         Axin(bcat=None, gsk3=None, ck1a=None, apc=ANY) + Gsk3(axin=None,dvl=None),
         kr_axin_gsk3)

    # APC unbinds
    Rule('Axin_unbinds_Apc',
         Axin(bcat=None, apc=1, gsk3=None, ck1a=None) % Apc(axin=1) >>
         Axin(bcat=None, apc=None, gsk3=None, ck1a=None) + Apc(axin=None),
         kr_axin_apc)

    Rule('Axin_Gsk3_unbinds_Apc',
         Axin(bcat=None, apc=1, gsk3=ANY, ck1a=None) % Apc(axin=1) >>
         Axin(bcat=None, apc=None, gsk3=ANY, ck1a=None) + Apc(axin=None),
         kr_axin_apc)

    Rule('Axin_Ck1a_unbinds_Apc',
         Axin(bcat=None, apc=1, gsk3=None, ck1a=ANY) % Apc(axin=1) >>
         Axin(bcat=None, apc=None, gsk3=None, ck1a=ANY) + Apc(axin=None),
         kr_axin_apc)

    # Bcat into dest comp
    # Here, we are requiring that the aa20 site of APC is NOT bound to another beta-catenin molecule in order for
    # beta-catenin to bind to Axin
    # NOTE: we are only allowing Bcat to bind to Axin if APC is unphosphorylated (this completes the loop and allows
    # the destruction complex to be recycled)
    Rule('Bcat_binds_DstCpx',
         Bcat(top=None, bottom=None, tcf4=None, loc='cyt', nterm='u') +
         Axin(bcat=None, ck1a=ANY, gsk3=2, apc=1) % Apc(axin=1, aa15=None, aa20=None, state='u') %
         Gsk3(axin=2, lithium=None, dvl=None) |
         Bcat(top=3, bottom=4, tcf4=None, loc='cyt', nterm='u') %
         Axin(bcat=3, ck1a=ANY, gsk3=2, apc=1) % Apc(axin=1, aa15=4, aa20=None, state='u') %
         Gsk3(axin=2, lithium=None, dvl=None),
         kf_bcat_dtcpx, kr_bcat_dtcpx)

    # We think Bcat can be phosphorylated and dephosphorylated when bound to Axin, whether it's bound to APC or not

    # NOTE: Bcat can only be bound to Axin at the 'top' site if ck1a, gsk3, and apc are all also bound
    # Therefore, we don't need to explicitly include those three sites in the rule
    # Also, dephosphorylation (reverse part of the rule) is implicitly modeled as due to PP2A
    Rule('Bcat_phos_Ck1a',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='u') % Axin(bcat=1) % Apc(aa15=2) |
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p1') % Axin(bcat=1) % Apc(aa15=2),
         k_bcat_phos_ck1a, k_dephos)  # ck1a phos first

    # GSK3 can't phosphorylate Bcat if it's bound to DLV or Li
    Rule('Bcat_phos_Gsk3',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p1') %
         Axin(bcat=1, gsk3=3) % Apc(aa15=2) % Gsk3(axin=3, lithium=None, dvl=None) >>
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p2') %
         Axin(bcat=1, gsk3=3) % Apc(aa15=2) % Gsk3(axin=3, lithium=None, dvl=None),
         k_bcat_phos_gsk3)  # then gsk3b phos

    # Don't restrict dephosphorylation of Bcat to case where Gsk3 is not bound to lithium
    # (i.e., Li could bind after phosphorylation of Bcat by Gsk3. Should be allowed to desphosphorylate in that case)
    Rule('Bcat_unphos_p2',
         Bcat(top=1, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1, gsk3=2) % Gsk3(axin=2, dvl=None) >>
         Bcat(top=1, tcf4=None, loc='cyt', nterm='p1') % Axin(bcat=1, gsk3=2) % Gsk3(axin=2, dvl=None),
         k_dephos)

    # # Bcat binds to aa15 site of APC (it is assumed phosphorylation state of Bcat does not affect binding)
    # Rule('Bcat_binds_aa15',
    #      Bcat(top=1, bottom=None, tcf4=None, loc='cyt') % Axin(bcat=1, apc=2) % Apc(axin=2, state='u', aa15=None) |
    #      Bcat(top=1, bottom=3, tcf4=None, loc='cyt') % Axin(bcat=1, apc=2) % Apc(axin=2, state='u', aa15=3),
    #      kf_bcat_apc, kr_bcat_apc)

    # APC phosphorylated by ck1a and gsk3
    # Dephosphorylation is assumed to be due to PP2A (implicit)
    Rule('Apc_phos_Ck1a',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1) %
         Apc(aa15=2, state='u', axin=ANY) % Ck1a() % Gsk3(dvl=None) >>
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1) %
         Apc(aa15=2, state='p1', axin=ANY) % Ck1a() % Gsk3(dvl=None),
         k_apc_phos_ck1a)

    Rule('Apc_unphos_p1', Apc(state='p1', aa20=None) >> Apc(state='u', aa20=None), k_dephos)

    Rule('Apc_phos_Gsk3',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1) %
         Apc(aa15=2, state='p1', axin=ANY) % Ck1a() % Gsk3(lithium=None, dvl=None) >>
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1) %
         Apc(aa15=2, state='p2', axin=ANY) % Ck1a() % Gsk3(lithium=None, dvl=None),
         k_apc_phos_gsk3)

    Rule('Apc_unphos_p2', Apc(state='p2', aa20=None) >> Apc(state='p1', aa20=None), k_dephos)

    # phosphorylation forces bcat to detach from axin
    Rule('Bcat_binds_Apc',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=1) %
         Apc(aa15=2, aa20=None, state='p2') % Ck1a() % Gsk3(dvl=None) >>
         Bcat(top=1, bottom=None, tcf4=None, loc='cyt', nterm='p2') % Axin(bcat=None) %
         Apc(aa15=None, aa20=1, state='p2') % Ck1a() % Gsk3(dvl=None),
         k_bcat_binds_apc_aa20)

    # apc is still bound to axin and bcat and apc are both phos

    #  Btrcp binds bcat
    # we are assuming that btrcp can bind (and unbind) to bcat whether it is ubiquitinated or not
    # Rule('btrcp_binds_bcat',
    #      Bcat(top=1, tcf4=None, loc='cyt', nterm='p2', bottom=None) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) +
    #      Btrcp(bcat=None) |
    #      Bcat(top=1, tcf4=None, loc='cyt', nterm='p2', bottom=2) % Apc(state='p2', aa20=1, aa15=None, axin=ANY) %
    #      Btrcp(bcat=2),
    #      kf_btrcp_binds_bcat, kr_btrcp_binds_bcat)

    # should axin=any be included?

    #  Btrcp ubiquitinates bcat
    # Rule('Bcat_ubiq',
    #      Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='x') % Apc(aa20=1) % Btrcp(bcat=2) >>
    #      Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='ub') % Apc(aa20=1) % Btrcp(bcat=2), k_bcat_ubiq)

    # Apc releases Bcat after ubiquitination
    Rule('Apc_releases_Bcat_ub',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt', state='ub') % Btrcp(b=2) % Apc(aa20=1) >>
         Bcat(top=None, bottom=2, tcf4=None, loc='cyt', state='ub') % Btrcp(b=2) + Apc(aa20=None),
         k_bcat_ub_release_apc)

    # Bcat degraded by proteosome
    # Rule('bcat_degradation', Bcat(top=None, bottom=2, tcf4=None, loc='cyt', state='ub') % Btrcp(bcat=2) >>
    # Btrcp(bcat=None), k_bcat_deg)

    Rule('Lithium_binds_GSK3',
         Gsk3(lithium=None, dvl=None) + Li(gsk3=None) | Gsk3(lithium=1, dvl=None) % Li(gsk3=1),
         kf_gsk3_li, kr_gsk3_li)

    return True


def create_wntmodel_rules(create=True):
    if not create:
        return False

    # # Beta Catenin binds btrcp
    # Rule('Bcat_binds_btrcp',
    #      Bcat(gsk3b=1, apc=2, btrcp=None, state='x', loc='cyt') % Gsk3b(bcat=1) % Apc(bcat=2) + Btrcp(b=None) | \
    #      Bcat(gsk3b=None, apc=None, btrcp=3, state='x', loc='cyt') % Btrcp(b=3) + Gsk3b(bcat=None) + Apc(bcat=None), \
    #      kf_btrcp_binds_bcat, kr_btrcp_binds_bcat)
    #
    # # Beta Catenin Ubiquitination
    # Rule('Bcat_Ubiq', Bcat(gsk3b=None, apc=None, btrcp=3, state='x', loc='cyt') % Btrcp(b=3) >>
    #      Bcat(gsk3b=None, apc=None, btrcp=3, state='ub', loc='cyt') % Btrcp(b=3), k_bcat_ubiq)
    #
    # # Beta catenin degradation
    # Rule('Bcat_degradation', Bcat(gsk3b=None, apc=None, btrcp=1, state='ub', loc='cyt') % Btrcp(b=1) >> Btrcp(b=None),
    #      k_bcat_deg)

    Rule('Bcat_binds_DVL',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt') % Axin(bcat=1, gsk3=3) % Gsk3(axin=3, dvl=None) % Apc(aa15=2) +
         Dvl(gsk3=None, rec=ANY) |
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt') % Axin(bcat=1, gsk3=3) % Gsk3(axin=3, dvl=4) % Apc(aa15=2) %
         Dvl(gsk3=4, rec=ANY),
         *k_bcat_dvl)

    # Release of APC from destruction complex by DVL
    Rule('DestCpx_Bcat_DVL_releases_APC',
         Bcat(top=1, bottom=2, tcf4=None, loc='cyt') %
         Axin(bcat=1, gsk3=3) % Gsk3(axin=3, dvl=4) % Dvl(gsk3=4, rec=ANY) % Apc(aa15=2) >>
         Bcat(top=1, bottom=None, tcf4=None, loc='cyt') %
         Axin(bcat=1, gsk3=3) % Gsk3(axin=3, dvl=4) % Dvl(gsk3=4, rec=ANY) + Apc(aa15=None),
         k_apc_release_dvl)

    # Beta catenin released from destruction complex and DVL
    # Assuming beta catenin dephosphorylates when it dissociates from the destruction complex
    Rule('DestCpx_DVL_releases_Bcat',
         Bcat(top=1, bottom=None, tcf4=None, nterm=WILD, loc='cyt') %
         Axin(bcat=1, gsk3=3) % Gsk3(axin=3, dvl=4) % Dvl(gsk3=4, rec=ANY) >>
         Bcat(top=None, bottom=None, tcf4=None, nterm='u', loc='cyt') +
         Axin(bcat=None, gsk3=None) + Gsk3(axin=None, dvl=None) + Dvl(gsk3=None, rec=ANY),
         k_bcat_release_dvl)

    # Beta catenin translocation to the nucleus
    Rule('Bcat_transport_nucleus',
         Bcat(top=None, bottom=None, loc='cyt') | Bcat(top=None, bottom=None, loc='nuc'),
         *k_bcat_nuc)

    # Monomer('Bcat', ['top', 'bottom', 'nterm', 'state', 'tcf4', 'loc'],  # Added tcf4 and loc components
    #         {'state': ['x', 'ub'], 'nterm': ['u', 'p1', 'p2'], 'loc': ['cyt', 'nuc']})

    # Beta catenin in nuclues binds to TCF4
    # TODO: Does tcf4 bind gGli2 before bcat binds?
    Rule('Bcat_binds_Tcf4',
         Bcat(tcf4=None, loc='nuc') + Tcf4(g=ANY, bcat=None) |
         Bcat(tcf4=1, loc='nuc') % Tcf4(g=ANY, bcat=1),
         *k_bcat_tcf4)

    # TCF4 binds to promoter of Gli2 gene
    Rule('Tcf4_binds_gGli2',
         Tcf4(g=None, bcat=None) + gGli2(tcf4=None) | Tcf4(g=1, bcat=None) % gGli2(tcf4=1),
         *k_tcf4_gli2)

    # SMAD3 binds to promoter of Gli2 gene
    Rule('Smad3_binds_gGli2',
         Smad3(g_gli2=None) + gGli2(smad3=None) | Smad3(g_gli2=1) % gGli2(smad3=1),
         *k_smad3_gli2)

    # Gli2 gene transcription
    Rule('gli2_transcription',
         gGli2(tcf4=1, smad3=ANY) % Tcf4(g=1, bcat=ANY) >>
         gGli2(tcf4=1, smad3=ANY) % Tcf4(g=1, bcat=ANY) + mGli2(),
         k_gli2_tx)

    # Gli2 translation
    Rule('gli2_translation',
         mGli2() >> mGli2() + Gli2(btrcp=None, g_pthlh=None, state='x', loc='cyt'),
         k_gli2_tl)

    # mGli2 degradation
    Rule('mGli2_degradation', mGli2() >> None, k_mgli2_deg)

    # TCF4 binds to promoter of PTHlH gene
    Rule('Tcf4_binds_gPthlh',
         Tcf4(g=None, bcat=None) + gPthlh(tcf4=None) | Tcf4(g=1, bcat=None) % gPthlh(tcf4=1),
         *k_tcf4_pthlh)

    # Gli2 binds to promoter of PTHlH gene
    Rule('Gli2_binds_gPthlh',
         Gli2(g_pthlh=None, loc='nuc') + gPthlh(gli2=None) |
         Gli2(g_pthlh=1, loc='nuc') % gPthlh(gli2=1),
         *k_gli2_pthlh)

    # PTHlH gene transcription
    Rule('gPthlh_transcription',
         gPthlh(tcf4=1, gli2=ANY) % Tcf4(g=1, bcat=ANY) >>
         gPthlh(tcf4=1, gli2=ANY) % Tcf4(g=1, bcat=ANY) + mPthlh(),
         k_pthlh_tx)

    # PTHrP translation
    Rule('Pthrp_translation', mPthlh() >> mPthlh() + Pthrp(), k_pthrp_tl)

    # mPthlh degradation
    Rule('mPthlh_degradation', mPthlh() >> None, k_mPthlh_deg)

    # PTHrP degrdation
    Rule('Pthrp_degradation', Pthrp() >> None, k_pthrp_deg)

    # Gli2 phosporylation
    Rule('Gli2_phosphorylation',
         Gli2(state='x', loc='cyt', btrcp=None) | Gli2(state='p', loc='cyt', btrcp=None),
         *k_gli2_phos)

    # Phopho-Gli2 translocates to nucleus
    Rule('Gli2_transport_nucleus',
         Gli2(state='p', g_pthlh=None, loc='cyt') | Gli2(state='p', g_pthlh=None, loc='nuc'),
         *k_gli2_nuc)

    # Gli2 binds BTRCP
    Rule('Gli2_binds_Btrcp',
         Gli2(state='x', loc='cyt', btrcp=None) + Btrcp(b=None) |
         Gli2(state='x', loc='cyt', btrcp=1) % Btrcp(b=1),
         *k_gli2_btrcp)

    # Gli2 ubiquitination
    Rule('Gli2_ubiquitination', Gli2(state='x', loc='cyt', btrcp=ANY) >> Gli2(state='ub', loc='cyt', btrcp=ANY),
         k_gli2_ubiq)

    # Gli2 degradation
    Rule('Gli2_degradation', Gli2(state='ub', btrcp=1) % Btrcp(b=1) >> Btrcp(b=None), k_gli2_deg)

    # WNT3A ligand binds receptor
    Rule('Wnt_binds_Rec', Wnt(rec=None) + Rec(wnt=None, dvl=None) | Wnt(rec=1) % Rec(wnt=1, dvl=None), *k_wnt_rec)

    # DVL binds WNT-bound receptor
    Rule('Dvl_binds_Wnt_Rec',
         Dvl(rec=None, gsk3=None) + Rec(wnt=1, dvl=None) % Wnt(rec=1) |
         Dvl(rec=2, gsk3=None) % Rec(wnt=1, dvl=2) % Wnt(rec=1),
         *k_dvl_rec)

    # DVL binds unbound receptor (only happens for high rigidity)
    # TODO: Commenting out the rule below for now, to simplify the debugging process. Will add it back in eventually
    # Rule('dvl_binds_rec', Dvl(rec=None) + Rec(wnt=None, dvl=None) | Dvl(rec=1) % Rec(wnt=None, dvl=1),
    # *k_dvl_rec_rigid)

    # DKK1 binds receptor
    Rule('Dkk1_binds_Rec', Dkk1(rec=None) + Rec(wnt=None) | Dkk1(rec=1) % Rec(wnt=1), *k_dkk1_rec)

    # WIF binds WNT3A
    Rule('Wif_binds_Wnt', Wif1(wnt=None) + Wnt(rec=None) | Wif1(wnt=1) % Wnt(rec=1), *k_wif_wnt)

    return True


create_destcpx_rules(create=True)
create_wntmodel_rules(create=True)

if __name__ == '__main__':

    # run simulation
    tspan = np.linspace(0, 1000, 1001)  #(0, 5000, 5001)  # (0, 40, 101)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)

    print('monomers %d' % len(model.monomers))
    print('rules %d' % len(model.rules))
    print('species %d' % len(model.species))
    print('reactions %d' % len(model.reactions))
    '''
    if True:
        plt.figure('gsk3_activity')
        Li_conc = np.arange(0, 101, 5)
        for k_phos in [1]:  # , 10, 100]:
            gsk3_activity = []
            for conc in Li_conc:
                print(conc, k_phos)
                result = sim.run(param_values={'Li_0': conc, 'k_bcat_phos_gsk3': k_phos})
                gsk3_activity.append(result.observables['GSK3_activity'][-1])
            gsk3_activity = np.array(gsk3_activity)
            plt.plot(Li_conc, gsk3_activity/gsk3_activity[0], 'o', label="k_phos=%g" % k_phos)
        plt.xlabel('Li concentration')
        plt.ylabel('GSK3 activity')
        plt.legend(loc=0)
    '''
    # run simulation(WNT)
    # if wntmodel_rules:
    if True:
        ncols = 2
        nrows = int(np.ceil(len(obs_to_plot) / ncols))
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, constrained_layout=True,
                                figsize=(6.4, 2.4*nrows))  # default: (6.4, 4.8)
        if len(axs.shape) == 1:
            axs = axs.reshape(axs.shape[0], 1)
        # run simulation
        traj = sim.run()
        # plot results
        row = 0
        col = 0
        for obs in obs_to_plot:
            print(obs)
            for o in obs:
                axs[col, row].plot(tspan, traj.observables[o], lw=2, label=o)
            axs[col, row].legend(loc=0)
            if (col+1) % ncols == 0:
                row += 1
                col = 0
            else:
                col += 1
        # add shared axis labels
        fig.supxlabel('Time')
        fig.supylabel('Concentration')
        # delete extra plots
        if col > 0:
            while col < ncols:
                fig.delaxes(axs[row][col])
                col += 1

    plt.show()
