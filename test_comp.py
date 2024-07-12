# # from mod import *
# # # from SynTemp.SynComp.valence_constrain import ChemicalRules
# # aldolAdd_F_id = ruleGMLString("""
# # rule [
# # 	#ruleID "Aldol Addition ->, id"
# # 	context [
# # 		edge [ source 1 target 2 label "=" ]
# # 		edge [ source 2 target 3 label "-" ]
# # 		edge [ source 3 target 4 label "-" ]
# # 		edge [ source 5 target 6 label "=" ]
# # 		node [ id 1 label "C" ]
# # 		node [ id 2 label "C" ]
# # 		node [ id 3 label "O" ]
# # 		node [ id 4 label "H" ]
# # 		node [ id 5 label "O" ]
# # 		node [ id 6 label "C" ]
# # 	]
# # ]
# # """)

# # aldolAdd_F = ruleGMLString("""
# # rule [
# # 	#ruleID "Aldol Addition ->"
# # 	left [
# # 		edge [ source 1 target 2 label "=" ]
# # 		edge [ source 2 target 3 label "-" ]
# # 		edge [ source 3 target 4 label "-" ]
# # 		edge [ source 5 target 6 label "=" ]
# # 	]
# # 	context [
# # 		node [ id 1 label "C" ]
# # 		node [ id 2 label "C" ]
# # 		node [ id 3 label "O" ]
# # 		node [ id 4 label "H" ]
# # 		node [ id 5 label "O" ]
# # 		node [ id 6 label "C" ]
# # 	]
# # 	right [
# # 		edge [ source 1 target 2 label "-" ]
# # 		edge [ source 2 target 3 label "=" ]
# # 		edge [ source 5 target 6 label "-" ]

# # 		edge [ source 4 target 5 label "-" ]
# # 		edge [ source 6 target 1 label "-" ]
# # 	]
# # ]
# # """)

# # r1 = aldolAdd_F_id
# # r2 = aldolAdd_F

# # p = GraphPrinter()
# # p.setReactionDefault()
# # p.withIndex = True
# # p.collapseHydrogens = False
# # r1.print(p)
# # r2.print(p)
# # from SynTemp.SynComp.rule_compose import RuleCompose
# # comp_rule = RuleCompose._compose(r1, r2)
# # print(len(comp_rule))
# # print(comp_rule[0].getGMLString())
# # def save_gml_from_text(gml_content, gml_file_path, rule_id, parent_ids):
# #     """
# #     Save a text string to a GML file by modifying the 'ruleID' line to include parent rule names.

# #     Parameters:
# #     - gml_content (str): The content to be saved to the GML file.
# #     - gml_file_path (str): The file path where the GML file should be saved.
# #     - rule_id (str): The original rule ID from the content.
# #     - parent_ids (list of str): List of parent rule IDs to prepend to the original rule ID.

# #     Returns:
# #     - bool: True if the file was successfully saved, False otherwise.
# #     """
# #     try:
# #         # Create the new ruleID by concatenating parent IDs with the original rule ID
# #         new_rule_id = 'p_' + '_'.join(parent_ids) + '_r_' + rule_id if parent_ids else 'r_' + rule_id

# #         # Initialize a list to hold the modified lines
# #         modified_lines = []

# #         # Iterate through each line and replace the 'ruleID' line as needed
# #         for line in gml_content.splitlines():
# #             if line.strip().startswith('ruleID'):
# #                 # Replace the whole line with the new ruleID
# #                 modified_lines.append(f'	ruleID "{new_rule_id}"')
# #             else:
# #                 modified_lines.append(line)

# #         # Join all lines back into a single string
# #         modified_content = "\n".join(modified_lines)

# #         # Write the modified content to the file
# #         with open(gml_file_path, "w") as file:
# #             file.write(modified_content)
# #         return True
# #     except FileNotFoundError:
# #         print(f"Unable to access the file path: {gml_file_path}")
# #         return False
# #     except Exception as e:
# #         print(f"An error occurred while writing to the file: {e}")
# #         return False

# # save_gml_from_text(comp_rule[0].getGMLString(), 'test_write.gml', rule_id = '0', parent_ids = ['0', '1'])
# # print("Computing all in MØD")
# # m = RCMatch(r1, r2)
# # modRes = m.composeAll()

# # def split(rs):
# # 	g = []
# # 	b = []
# # 	for r in rs:
# # 		if checkRule(r):
# # 			g.append(r)
# # 		else:
# # 			b.append(r)
# # 	return g, b

# # valence_check = ChemicalRules()
# # goodMod, badMod = valence_check.split(modRes)

# # print("|MØD|:", len(goodMod) + len(badMod))
# # print("  |good|:", len(goodMod))
# # print("  |bad|:", len(badMod))

# # postSection("Good")
# # for r in goodMod:
# # 	r.print()
# # postSection("Bad")
# # for r in badMod:
# # 	r.print()
# from SynTemp.SynComp.rule_compose import RuleCompose
# #single_rule_path = 'SingleRule'
# single_rule_path = 'Data/DPO/USPTO_50K/Hydrogen/R0'
# compose_path = 'Compose'

# #RuleCompose._auto_compose(rule_path=single_rule_path, rule_path_compose=compose_path)