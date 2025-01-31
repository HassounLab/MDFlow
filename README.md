# Analysis of Metabolic Network Disruption in Engineered Microbial Hosts due to Enzyme Promiscuity

* Vladimir Porokhin, Department of Computer Science, Tufts University, Medford, MA, vladimir.porokhin@tufts.edu
* Sara A. Amin, Department of Computer Science, Tufts University, Medford, MA, sara.amin@tufts.edu
* Trevor B. Nicks, Department of Chemical and Biological Engineering, Tufts University, Medford, MA, trevor.nicks@tufts.edu
* Venkatesh Endalur Gopinarayanan, Department of Chemical and Biological Engineering, Tufts University, Medford, MA, venkatesh.endalur_gopinarayanan@tufts.edu
* Nikhil U. Nair, Department of Chemical and Biological Engineering, Tufts University, Medford, MA, nikhil.nair@tufts.edu
* Soha Hassoun, Department of Computer Science and Department of Chemical & Biological Engineering, Tufts University, Medford, MA, soha.hassoun@tufts.edu

# Multicopy Suppressors

Multicopy suppressor cases are analyzed using the `patrick-fba.py` script. It requires the following files to be present in the `data` directory:

* `reactions-patrick.tsv` (included): list of promiscuous reactions generated by PROXIMAL
* `iML1428-iso_Glucose.json`: iML1428 model from Supplementary Dataset 1 in [1]
* `kegg-list-eco.txt`: list of *E. coli* genes from [KEGG](http://rest.kegg.jp/list/eco).
* `kegg-link-eco-ec.txt`: list of *E. coli* enzymes linked to genes, also from [KEGG](http://rest.kegg.jp/link/eco/ec).

The script should be run twice, using FBA and MOMA separately:

	python patrick-fba.py fba
	python patrick-fba.py moma

These commands can be executed in parallel. The script generates tables as in the Supplementary File 1 accompanying the paper and writes them to the `output` directory.

[1]: Monk, J., Lloyd, C., Brunk, E. *et al.* *i*ML1515, a knowledgebase that computes *Escherichia coli* traits. *Nat Biotechnol* **35**, 904-908 (2017). https://doi.org/10.1038/nbt.3956

# 3-HP Synthesis Pathways

Processing the 3-HP synthesis pathways is done using two scripts. The first one, `3hp-fba.py`, is responsible for conducting the simulations of the baseline, engineered, and disrupted models following the method described in the paper. The following files must to be present in the `data` directory to run the script:

* `data/reactions-3hp.tsv` (included): list of promiscuous reactions generated by PROXIMAL
* `data/iML1515.json`: iML1515 model from the [BiGG website](http://bigg.ucsd.edu/models/iML1515).

The script needs to be run three times for each pathway (`p1` and `p2`): for scenario 1 (`s1`), scenario 2 (`s2`), and both combined (`fixed`). Thus, a full analysis requires running six commands as follows:

	python 3hp-fba.py p1 s1 fba %seed%
	python 3hp-fba.py p1 s2 fba %seed%
	python 3hp-fba.py p1 fixed fba %seed%
	python 3hp-fba.py p2 s1 fba %seed%
	python 3hp-fba.py p2 s2 fba %seed%
	python 3hp-fba.py p2 fixed fba %seed%

These commands can all be executed in parallel to maximize CPU utilization. Results of the simulations will be placed in the `output` directory.

Once all simulations are completed, `3hp-fig.py`, can be used to compile the results into scatter plots similar to the ones presented in the paper:

	python 3hp-fig.py

The figures will be placed in the `output` directory.
