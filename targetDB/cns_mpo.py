#!/usr/bin/env python
# ----------------------------------------------------------------------#
# 			CNS MPO Score and Solubility forecaster index (3.2)			#
# --------------------------- Compatibility ----------------------------#
# 					Compatible with Python 2.7 & 3.5					#
# ----------------------Prerequiste installations-----------------------#
# 																		#
# 		- Chemaxon Marvin suite WITH license							#
# 				cxcalc module required to work 	(used for pKa pred)		#
# 					(make sure to add cxcalc to $PATH)					#
# by default : C:\Program Files (x86)\ChemAxon\MarvinBeans\bin			#
# How to use :															#
# 		- input the sdf file name when requested (.sdf included)		#
# 		- the script output an sdf with #name_out.sdf					#
# 			This sdf contains the fields with : 						#
# 						- CNS MPO Score									#
# 						- Solubility forecaster index (SFI)[optional]	#
# 						- bpKa,logD(7.4),logP,MW,HBD,#Ar,TPSA			#
# 						- All individual components of MPO Score		#
# ----------------------------------------------------------------------#
# (C) Dr. De Cesco Stephane - v3.2 - 22/03/2017							#
# ----------------------------------------------------------------------#
# -------------------------- MODULES IMPORT ----------------------------#

def monotonic_score(value, lower, upper):
    # =======================================================================================
    #         lower
    # 1|---------
    #  |		  \
    #  |           \@value
    #  |			\
    # 0|_____________\______________________
    # 	             upper
    # Function to return a score between 0 and 1 depending of the value of the parameters.
    # 		 | 1 if value < lower
    # Score ={ f(upper and lower) if upper < value < lower
    # 		 | 0 if value > upper
    # =======================================================================================
    upper = float(upper)
    lower = float(lower)
    v = value.copy()
    v[value<= lower] = 1
    v[value >= upper] = 0
    v[(value > lower) & (value < upper)] = 1 - ((value - lower) * (1 / (upper - lower)))

    return v


def hump_score(value, low1, up1, up2, low2):
    # =======================================================================================
    #         	up1		 up2
    # 1|		  --------
    #  |		 / 		  \
    #  |        /      	   \
    #  |	   /			\
    # 0|______/______________\_______________
    # 	     low1			  low2
    # Function to return a score between 0 and 1 depending of the value of the parameters.
    # 		 | 0 if value < low1
    # 		 | f (low1 and up1) if low1 < value < up1
    # Score ={ 1 if up1 < value < up2
    # 		 | f (up2 and low2) if up2 < value < low2
    # 		 | 0 if value > lower
    # =======================================================================================

    low1, up1, up2, low2 = float(low1), float(up1), float(up2), float(low2)
    v = value.copy()
    v[value <= low1] = 0
    v[value > low2] = 0
    v[(up1 < value)&(value<= up2)] = 1
    v[(low1 < value) & (value<= up1)]= ((value - low1) * (1 / (up1 - low1)))
    v[(up2 < value) & (value<= low2)]= 1 - ((value - up2) * (1 / (low2 - up2)))
    return v


def calc_mpo_score(bpka=None,logP=None,logD=None,MW=None,HBD=None,TPSA=None):
    bpKa_score = monotonic_score(bpka, 8,10)
    logP_score = monotonic_score(logP, 3, 5)
    logD_score = monotonic_score(logD, 2, 4)
    MW_score = monotonic_score(MW, 360, 500)
    HBD_score = monotonic_score(HBD, 0.5, 3.5)
    TPSA_score = hump_score(TPSA, 20, 40, 90, 120)
    return bpKa_score + logP_score + logD_score + MW_score + HBD_score + TPSA_score

