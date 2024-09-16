from mod import *

rule1 = """
rule [
   ruleID "8"
   left [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "N" ]
      node [ id 2 label "C" ]
      node [ id 3 label "H" ]
      node [ id 4 label "O" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
]
"""

rule2 = """
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

rule3 = """
rule [
   ruleID "6"
   left [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "N" ]
      node [ id 2 label "H" ]
      node [ id 3 label "C" ]
      node [ id 4 label "O" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
]
"""
rule4 = """
rule [
   ruleID "5"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "O" ]
      node [ id 3 label "H" ]
      node [ id 4 label "O" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]
"""

rule5 = """
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
rule6 = """
rule [
   ruleID "12"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "B" ]
      node [ id 2 label "C" ]
      node [ id 3 label "C" ]
      node [ id 4 label "Br" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]
"""

rule7 = """
rule [
   ruleID "14"
   left [
      edge [ source 1 target 2 label "=" ]
      edge [ source 3 target 6 label "-" ]
      edge [ source 4 target 5 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "O" ]
      node [ id 3 label "N" ]
      node [ id 4 label "H" ]
      node [ id 5 label "H" ]
      node [ id 6 label "H" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 5 label "-" ]
      edge [ source 2 target 6 label "-" ]
   ]
]
"""

rule8 = """
rule [
   ruleID "3"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "H" ]
      node [ id 2 label "O" ]
      node [ id 3 label "C" ]
      node [ id 4 label "Cl" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]
"""

rule9 = """
rule [
   ruleID "1"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "Br" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]
"""
rule10 = """
rule [
   ruleID "13"
   left [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "S" ]
      node [ id 2 label "H" ]
      node [ id 3 label "Cl" ]
      node [ id 4 label "N" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]
"""
rule11 = """
rule [
   ruleID "10"
   left [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
   context [
      node [ id 1 label "I" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
]
"""
rule12 = """
rule [
   ruleID "21"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "I" ]
      node [ id 3 label "H" ]
      node [ id 4 label "N" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]
"""
rule13 = """
rule [
   ruleID "4"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "=" ]
   ]
   context [
      node [ id 1 label "H" ]
      node [ id 2 label "N" ]
      node [ id 3 label "N" ]
      node [ id 4 label "C" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
]
"""
rule14 = """
rule [
   ruleID "23"
   left [
      edge [ source 1 target 2 label "=" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "O" ]
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
rule15 = """
rule [
   ruleID "17"
   left [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "N" ]
      node [ id 3 label "H" ]
      node [ id 4 label "F" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
]
"""
rule16 = """
rule [
   ruleID "7"
   left [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "=" ]
   ]
   context [
      node [ id 1 label "H" ]
      node [ id 2 label "C" ]
      node [ id 3 label "C" ]
      node [ id 4 label "H" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""

rule17 = """
rule [
   ruleID "55"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      edge [ source 4 target 5 label "-" ]
      edge [ source 6 target 7 label "-" ]
   ]
   context [
      node [ id 1 label "O" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
      node [ id 5 label "H" ]
      node [ id 6 label "H" ]
      node [ id 7 label "H" ]
   ]
   right [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 2 target 6 label "-" ]
      edge [ source 2 target 7 label "-" ]
      edge [ source 3 target 5 label "-" ]
   ]
]
"""

rule18 = """
rule [
   ruleID "44"
   left [
      edge [ source 1 target 4 label "-" ]
      edge [ source 1 target 5 label "-" ]
      edge [ source 2 target 3 label "=" ]
   ]
   context [
      node [ id 1 label "N" ]
      node [ id 2 label "O" ]
      node [ id 3 label "C" ]
      node [ id 4 label "H" ]
      node [ id 5 label "H" ]
   ]
   right [
      edge [ source 1 target 3 label "=" ]
      edge [ source 2 target 4 label "-" ]
      edge [ source 2 target 5 label "-" ]
   ]
]
"""
rule19 = """
rule [
   ruleID "46"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]
   context [
      node [ id 1 label "Cl" ]
      node [ id 2 label "C" ]
      node [ id 3 label "B" ]
      node [ id 4 label "C" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""
rule20 = """
rule [
   ruleID "15"
   left [
      edge [ source 1 target 4 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
   context [
      node [ id 1 label "H" ]
      node [ id 2 label "S" ]
      node [ id 3 label "Cl" ]
      node [ id 4 label "O" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
   ]
]
"""

for i in range(1, 21):
    rule_var = locals()[f"rule{i}"]  # Access the local variable dynamically
    ruleGMLString(rule_var)  # Call the function with the dynamically accessed variable


for rule in inputRules:
    rule.print(second=True)
