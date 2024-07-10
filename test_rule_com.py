from SynTemp.SynUtils.graph_utils import load_gml_as_text
from mod import *

rule_0 = load_gml_as_text('Data/DPO/USPTO_50K/Hydrogen/R0/0.gml')
rule_1 = load_gml_as_text('Data/DPO/USPTO_50K/Hydrogen/R0/1.gml')

reaction_rule_0 = ruleGMLString(rule_0, invert=False, add=True)
reaction_rule_1 = ruleGMLString(rule_1, invert=False, add=True)


rc = rcEvaluator(inputRules)
# The special global object 'rcParallel' is used to make a pseudo-operator:
exp = rcId(reaction_rule_0) *rcParallel*  rcUnbind(reaction_rule_1)
rules = rc.eval(exp)
for p in rules:
   p.print()