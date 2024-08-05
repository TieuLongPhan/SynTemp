from mod import *

rule0 = """
rule [
   ruleID "Reduction Alkyne to Alkene"
   left [
      edge [ source 1 target 2 label "#" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "H" ]
      node [ id 4 label "H" ]
   ]
   right [
      edge [ source 1 target 2 label "=" ]
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""

rule1 = """
rule [
   ruleID "Reduction Alkene to Alkane"
   left [
      edge [ source 1 target 2 label "=" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "H" ]
      node [ id 4 label "H" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""

rule01 = """
rule [
	ruleID "Reduction Alkyne to Alkane"
	labelType "string"
	left [
		edge [ source 0 target 1 label "#" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 4 target 5 label "-" ]
	]
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		node [ id 2 label "H" ]
		node [ id 3 label "H" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
	]
	right [
		edge [ source 0 target 1 label "-" ]
		edge [ source 0 target 2 label "-" ]
		edge [ source 0 target 4 label "-" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 5 label "-" ]
	]
]
"""


for rule in [rule0, rule1, rule01]:

    ruleGMLString(rule)


for rule in inputRules:
    rule.print()
