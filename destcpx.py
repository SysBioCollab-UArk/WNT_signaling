from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator

Model()

#monomer
Monomer('Axin', ['bcat', 'gsk3', 'ck1a', 'apc'])
Monomer('GSK3', ['axin'])
Monomer('CK1A', ['axin'])
Monomer('APC', ['axin', 'aa15', 'aa20'],{'aa20': ['u', 'p']})
Monomer('Bcat', ['top', 'bottom', 'nterm', 'btrcp'], {'btrcp': ['x', 'ub'],'nterm':['u', 'p1', 'p2']}) # p1 by ck1a , p2 by gsk3
Monomer('BTRCP', ['bcat'])

#Initials
Parameter('Axin_0', 10)
Parameter('GSK3_0', 50)
Parameter('CK1A_0', 50)
Parameter('APC_0', 50)
Parameter('Bcat_0', 100)
Parameter('BTRCP_0', 50)

Initial(Axin(bcat=None, gsk3=None, ck1a=None, apc=None), Axin_0)
Initial(GSK3(axin=None), GSK3_0)
Initial(CK1A(axin=None), CK1A_0)
Initial(APC(axin=None, aa15=None, aa20='u'), APC_0)
Initial(Bcat(top=None, bottom=None, nterm='u', btrcp='x'), Bcat_0)
Initial(BTRCP(bcat=None), BTRCP_0)

#observables
Observable('Axin_tot', Axin())#total amount of axin
Observable('GSK3_tot', GSK3())#total amount of gsk3
Observable('CK1A_tot', CK1A())#total amount of ck1a
Observable('APC_aa20u', APC(aa20='u'))
Observable('APC_aa20p', APC(aa20='p'))
Observable('Bcat_u', Bcat(nterm='u'))
Observable('Bcat_p1', Bcat(nterm='p1'))
Observable('Bcat_p2', Bcat(nterm='p2'))
Observable('Bcat_ub', Bcat(btrcp='ub'))
Observable('BTRCP_tot', BTRCP())


#Rules

#Axin binding rules
Parameter('kf_axin_ck1a',1)
Parameter('kf_axin_gsk3',1)
Parameter('kf_axin_apc',1)
Rule('axin_binds_ck1a', Axin(ck1a=None) + CK1A(axin=None) >> Axin(ck1a=1) % CK1A(axin=1), kf_axin_ck1a)
Rule('axin_binds_gsk3', Axin(gsk3=None) + GSK3(axin=None) >> Axin(gsk3=1) % GSK3(axin=1), kf_axin_gsk3)
Rule('axin_binds_apc', Axin(apc=None) + APC(axin=None) >> Axin(apc=1) % APC(axin=1), kf_axin_apc)

#Bcat entering the destruction complex
Parameter('kf_bcat_dtcpx', 1)
Parameter('kf_bcat_apc', 1)
Rule('Bcat_binds_dtcpx', Bcat(top=None) + Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=ANY) >> \
     Bcat(top=1) % Axin(bcat=1, ck1a=ANY, gsk3=ANY, apc=ANY), kf_bcat_dtcpx )
Rule('Bcat_binds_aa15', Bcat(top=ANY,bottom=None) + APC(aa15=None) >> Bcat(top=ANY, bottom=1) % APC(aa15=1), kf_bcat_apc)
#Bcat phosphorilated by gsk3 and ck1a
Parameter('kf_bcat_phos_gsk3', 1)
Parameter('kr_bcat_u_gsk3', 1)
Parameter('kf_bcat_phos_ck1a', 1)
Parameter('kf_bcat_u_ck1a', 1)
Rule('bcat_p_gsk3', Bcat(nterm='u', top=ANY) | Bcat(nterm='p1', top=ANY), kf_bcat_phos_gsk3, kr_bcat_u_gsk3)
Rule('bcat_p_ck1a', Bcat(nterm='u', top=ANY) | Bcat(nterm='p2', top=ANY), kf_bcat_phos_ck1a, kr_bcat_u_ck1a)
#APC phosphorilated by ck1a
Parameter('kf_phos_ck1a', 2)
Parameter('kr_phos_ck1a', 2)
Rule('apc_p_ck1a', APC(aa20='u') | APC(aa20='p'), kf_phos_ck1a,kr_phos_ck1a)
#APC binds to bcat
Parameter('kf_bcat_binds_apc', 1)
Parameter('kr_bcat_binds_apc', 1)
Rule('bcat_binds_apc', Bcat(top=1) % Axin(bcat=1, ck1a=ANY, gsk3=ANY, apc=ANY) + APC(aa20=None) | \
     Bcat(top=1) + Axin(bcat=None, ck1a=ANY, gsk3=ANY, apc=ANY) % APC(aa20=1), kf_bcat_binds_apc, kr_bcat_binds_apc)

#testing commit
