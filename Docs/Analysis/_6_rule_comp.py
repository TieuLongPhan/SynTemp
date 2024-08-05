from mod import *

p_42 = """
rule [
   ruleID "42"
   left [
      edge [ source 1 target 4 label "=" ]
      edge [ source 2 target 3 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "H" ]
      node [ id 4 label "O" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 1 target 4 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
]
"""

p_94 = """
rule [
   ruleID "94"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "O" ]
      node [ id 2 label "H" ]
      node [ id 3 label "C" ]
      node [ id 4 label "S" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""

p_42_94_r_0 = """
rule [
   ruleID "163"
   left [
      edge [ source 1 target 2 label "=" ]
      edge [ source 3 target 4 label "-" ]
      edge [ source 3 target 5 label "-" ]
   ]
   context [
      node [ id 1 label "O" ]
      node [ id 2 label "C" ]
      node [ id 3 label "C" ]
      node [ id 4 label "S" ]
      node [ id 5 label "H" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 4 target 5 label "-" ]
   ]
]
"""


p_2 = """
rule [
   ruleID "2"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "Cl" ]
      node [ id 3 label "N" ]
      node [ id 4 label "H" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""

p_0 = """
rule [
   ruleID "0"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "N" ]
      node [ id 2 label "H" ]
      node [ id 3 label "C" ]
      node [ id 4 label "Br" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""

p_0_2_r_1 = """
rule [
   ruleID "227"
   left [
      edge [ source 1 target 6 label "-" ]
      edge [ source 1 target 7 label "-" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 4 target 5 label "-" ]
   ]
   context [
      node [ id 1 label "N" ]
      node [ id 2 label "C" ]
      node [ id 3 label "Cl" ]
      node [ id 4 label "C" ]
      node [ id 5 label "Br" ]
      node [ id 6 label "H" ]
      node [ id 7 label "H" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 7 label "-" ]
      edge [ source 5 target 6 label "-" ]
   ]
]
"""

for rule_var in [p_42, p_94, p_42_94_r_0, p_0, p_2, p_0_2_r_1]:
    ruleGMLString(rule_var)


for rule in inputRules:
    rule.print(second=True)
