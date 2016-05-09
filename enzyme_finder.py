import PyFBA
import sys

compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes('gramnegative')

looking_for = 'inosine'
def reaction_search(looking_for):
	for c in compounds:
		if looking_for in c.lower():
			print("Found " + c + " Reactions:  " + str(compounds[c].reactions))

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
     	
     	
""" 
reactions['rxn00516'].enzymes
set(['cpx00362', 'cpx34550', 'cpx.357'])
>>>enzymes["cpx.362"].roles
set(['Short chain fatty acids transporter'])

Found 7-methylinosine (location: ) Reactions:  set([])
Found Inosine 5'-tetraphosphate (location: ) Reactions:  set([])
Found tRNA-Arg-inosine34 (location: ) Reactions:  set([])
Found Lentiginosine (location: ) Reactions:  set([])
Found 5'-Oxoinosine (location: c) Reactions:  set(['rxn00042', 'rxn11662', 'rxn04055'])
Found Inosine (location: e) Reactions:  set(['rxn08754', 'rxn29210', 'rxn28604', 'rxn05316', 'rxn08753'])
Found 5'-Butyrylphosphoinosine (location: c) Reactions:  set(['rxn01238'])
Found 6-Thioinosine 5'-monophosphate (location: ) Reactions:  set([])
Found 6-Thioinosine-5'-diphosphate (location: c) Reactions:  set(['rxn11770', 'rxn12025', 'rxn11769', 'rxn12026'])
Found tRNA-Arg-inosine34 (location: c) Reactions:  set(['rxn37808', 'rxn25976'])
Found Inosine (location: c) Reactions:  set(['rxn37600', 'rxn01137', 'rxn19933', 'rxn30598', 'rxn00835', 'rxn30713', 'rxn12445', 'rxn01358', 'rxn27951', 'rxn37599', 'rxn29210', 'rxn19932', 'rxn00042', 'rxn27949', 'rxn30907', 'rxn29985', 'rxn32620', 'rxn15248', 'rxn20983', 'rxn31938', 'rxn32622', 'rxn05316', 'rxn31940', 'rxn04054', 'rxn32633', 'rxn31942', 'rxn08754', 'rxn00831', 'rxn31888', 'rxn01299', 'rxn28604', 'rxn18974', 'rxn08753'])
Found 6-Thioinosine 5'-monophosphate (location: c) Reactions:  set(['rxn12025', 'rxn11771', 'rxn11767', 'rxn11765', 'rxn11769', 'rxn11768'])
Found 1-methylinosine (location: ) Reactions:  set([])
Found Inosine (location: h) Reactions:  set(['rxn12445'])
Found Inosine-5'-carboxylate (location: ) Reactions:  set([])
Found 5'-S-Methyl-5'-thioinosine (location: ) Reactions:  set([])
Found Inosine (location: ) Reactions:  set([])
Found 5'-Oxoinosine (location: ) Reactions:  set([])
Found 6-Thioinosine-5'-triphosphate (location: c) Reactions:  set(['rxn11770', 'rxn11771', 'rxn12026'])
Found Deoxyinosine (location: c) Reactions:  set(['rxn01858', 'rxn30004', 'rxn09010', 'rxn08343', 'rxn30042', 'rxn11225', 'rxn01985', 'rxn29209', 'rxn08344', 'rxn32501', 'rxn19457', 'rxn10055', 'rxn32504', 'rxn30597', 'rxn32204', 'rxn32203', 'rxn09687', 'rxn12417'])
Found 6-Methylthioinosine-5'-monophosphate (location: c) Reactions:  set(['rxn11767', 'rxn11766'])
Found 5'-Acylphosphoinosine (location: ) Reactions:  set([])
Found Deoxyinosine (location: ) Reactions:  set([])
Found 6-Thioinosine-5'-diphosphate (location: ) Reactions:  set([])
Found Inosine 5'-tetraphosphate (location: c) Reactions:  set(['rxn00516'])
Found 5'-Butyrylphosphoinosine (location: ) Reactions:  set([])
Found 5'-S-Methyl-5'-thioinosine (location: c) Reactions:  set(['rxn16511', 'rxn16503'])
Found Deoxyinosine (location: e) Reactions:  set(['rxn08343', 'rxn29209', 'rxn09687', 'rxn08344'])
Found 6-Methylthioinosine-5'-monophosphate (location: ) Reactions:  set([])
Found Deoxyinosine (location: h) Reactions:  set(['rxn12417'])
Found 6-Thioinosine-5'-triphosphate (location: ) Reactions:  set([])
"""