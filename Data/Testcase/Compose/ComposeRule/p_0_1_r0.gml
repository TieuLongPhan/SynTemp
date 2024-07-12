rule [
	ruleID "p_0_1_r_0"
	labelType "string"
	left [
		edge [ source 0 target 1 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 6 target 7 label "-" ]
	]
	context [
		node [ id 0 label "H" ]
		node [ id 1 label "N" ]
		node [ id 2 label "C" ]
		node [ id 3 label "Br" ]
		node [ id 4 label "Br" ]
		node [ id 5 label "C" ]
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
	]
	right [
		edge [ source 0 target 3 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 4 target 7 label "-" ]
		edge [ source 5 target 6 label "-" ]
	]
]