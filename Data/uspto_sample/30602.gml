rule [
   ruleID "30602"
   left [
      edge [ source 2 target 3 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "H" ]
      node [ id 3 label "N" ]
   ]
   right [
      edge [ source 1 target 3 label "-" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "-" ]
   ]
]