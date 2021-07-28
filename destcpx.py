from pysb import *
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator

Model()

#monomer
Monomer('Axin', ['Bcat', 'gsk3', 'ck1a', 'apc'])
Monomer('GSK3', ['axin', 'apc','bcat'])
Monomer('CK1A', ['axin', 'apc', 'bcat')
Monomer('APC', ['axin', 'bcat', 'ck1a', 'gsk3'], {'ck1a':['p', 'np'], 'gsk3': ['p', 'np']})#ASK IF THE STATE OF THESE MONOMERS ARE NECESSARY
Monomer('Bcat', ['axin', 'apc', 'ck1a', 'gsk3', 'btrcp'], {'ck1a':['p', 'np'], 'gsk3': ['p', 'np'], 'apc': ['aa15', 'aa20'], 'btrcp': ['x', 'ub']})# CK1A site can be phosphorylated or not phosphorylated as well as GSK3 site
Monomer('BTRCP', ['bcat'])

#Initials
Parameter('Axin_0', 10)
Parameter('GSK3_0', 50)
Parameter('CK1A_', 50)
Parameter('APC_0', 50)
Parameter('Bcat_0', 100)
Parameter('BTRCP', 50)

Initial(Axin(bcat=None, gsk3=None, ck1a=None, apc=None), Axin_0)
Initial(GSK3(axin=None, apc=None,bcat=None), GSK3_0)
Initial(CK1A(axin=None, apc=None, bcat=None), CK1A_0)
Initial(APC(axin=None, bcat=None, ck1a='np', gsk3='np') % Bcat(apc=1), APC_0)
Initial(Bcat(axin=None, apc=None, ck1a='np', gsk3='np', btrcp='x'), Bcat_0)
Initial(BTRCP(bcat=None))

#observables
Observable('Axin_tot', Axin())#total amount of axin
Observable('GSK3_tot', Gsk3())#total amount of gsk3
Observable('CK1A_tot', CK1A())#total amount of ck1a
Observable('APC_PHO_BY_CK1A', APC(ck1a='p'))# APC phosphorilated by ck1a
Observable('APC_PHO_BY_GSK3', APC(gsk3='p'))# APC phosphorilated by gsk3
Observable('Bcat_PHO_BY_CK1A', Bcat(ck1a='p'))# Bcat phosphorilated by ck1a
Observable('Bcat_PHO_BY_GSK3', Bcat(gsk3='p'))# Bcat phosphorilated by gsk3
Observable('Bcat_APC_boundA', Bcat(apc='aa15'))#Bcat bound to APC at replicates 15aa
Observable('Bcat_APC_boundA', Bcat(btrcp='ub'))#Bcat ubiquitinated by btrcp
Observable('BTRCP_tot', BTRCP())

#Rate constants
Parameter('k_gsk3_ckia_axin', 100)# rate at which Axin binds to gsk3 and ck1a
Parameter('Bcat_join_the_complex', 150)

#Rules

#Axin binds to gsk3 and ck1a
Rule('Axin_bindind', Axin(gsk3b=None,ck1a=None, apc=None) % GSK3B(axin=None) % CK1A(axin=none) % APC(axin=None) >> \
     Axin(gsk3b=3, ck1a=6,apc=12) % GSK3(axin=3) % CK1A(axin=6),APC(axin=12), k_gsk3_ckia_axin)# not sure if it can be reversible

#Bcat entering the destruction complex
Rule('Bcat_Dcomplex',Bcat(axin= None, apc=None) % Axin(bcat=None, gsk3b=3, ck1a=6) % APC(bcat=None,axin=12) >> \
     Bcat(axin=7,apc='15aa') % Axin(bcat=7, gsk3b=3, ck1a=6) % APC(bcat=7, axin=12), Bcat_join_the_complex )# still confused about the site at where APC(binds Bcat)
#GK1A phosphorilates Bcat
Rule('Bcat_phos_ck1a', Bcat(ck1a='np',axin='4',apc='15aa') % CK1A(axin=9,bcat=None) )