# Written by Blaire Robinson
# 17 March 2016
#
# Description: uses PyFBA database to search for enzymes associated with compounds of interest.

import PyFBA
import sys

compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes('gramnegative')

# variable looking_for can be changed for compound of interest
looking_for = 'inosine'

# fuction that searches of all reactions associated with specified compound
def reaction_search(looking_for):
	for c in compounds:
		if looking_for in c.lower():
			print("Found " + c + " Reactions:  " + str(compounds[c].reactions))

# function that searches for all enzymes associated with specified compound.
def enzyme_search(c):
	for c in compounds:	
		if looking_for in c.lower():
			rxns_set = compounds[c].reactions
			for x in rxns_set:
				enzyme_set = reactions[x].enzymes
				for e in enzyme_set:
					try:
						print("Enzyme: " + e + "; Role: " + str(enzymes[e].roles))
					except KeyError:
						print >> sys.stderr, e + "\t" + x
     	
reaction_search(looking_for)
enzyme_search(looking_for)
