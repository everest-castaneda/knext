
KEGG NetworkX Topological (KNeXT) parser
========================================

KNeXT downloads and parses Kyoto Encylopedia of Genes and Genomes 
(KEGG) markup language files (KGML). This tool employs NetworkX's framework
to create gene-only networks, but mixed (gene, compound, pathway) networks
can also be generated. All output files are in TSV format. KNeXT also
retrieves a TXT file of node x-y axis coordinates for use in NetworkX's
graph visualization library, and it is able to convert KEGG IDs 
into Uniprot and NCBI IDs. KNeXT also maximizes metadata information
through preserving each edge's information and adding gene, compound, and
pathways names.

Usage
-----

.. code:: text

    Primary line: knext get-kgml [SPECIES_NAME]
      
      KEGG NetworkX Topological (KNeXT) parser uses the KEGG
      API to gather all KGML files for a single species. 
      Input species name in 3 to 4 letter KEGG organism code. 
    
    Options:
      --help,	shows options and website for KEGG organism codes
      -d/--d,	directory in which to save output

    Primary line: knext genes [Input]

      KNeXT parser deploy's NetworkX's
      framework to create gene-only representations of KGML files.
      Genes between compounds are propagated before compounds are dropped.

    Options:
      Input	KGML file or folder of KGML files to parse
      -r/--results	file or folder where output should be stored	
      -g/--graphics	outputs TXT file of x-y axis coordinates
      -u/--unique	TSV file's genes have a terminal modifier
      -n/--names    TSV file includes a column with gene names
		    Notice: increases parsing time and internet connection required
      --help	shows options and file types

    Primary line: knext mixed [Input]

      KNeXT parser creates mixed (genes, compounds, pathways)
      representations of KGML files.

    Options:
      Input	KGML file or folder of KGML files to parse
      -r/--results	file or folder where output should be stored
      -g/--graphics	outputs TXT file of x-y axis coordinates
      -u/--unique	TSV file's genes have a terminal modifier
      -n/--names TSV file includes a column with compound, gene, and pathway names
		 Notice: increases parsing time and internet connection required
      --help	shows options and file types

    Primary line: knext convert [OPTIONS]
      
      KNeXT parser converts KEGG entry IDs in TSV output files into
      UniProt or NCBI IDs.
    
    Options:
      file	PATH:	path to TSV file
      species	TEXT:	KEGG 3 to 4 letter organism code
      --uniprot	optional flag for output:	use if UniProt IDs are the desired output
      --unique	optional flag for output:	use if the TSV file has terminal modifiers
      --graphics	PATH:	graphics file
      --help	optional flag:	shows options

    Options:
      folder	PATH:	path to folder containing TSV files         
      species	TEXT:	KEGG 3 to 4 letter organism code
      --uniprot	optional flag for output:         use if UniProt IDs are the desired output
      --unique	optional flag for output:         use if the TSV file has terminal modifiers   
      --graphics	PATH:       path to folder containing graphics files          
      --help	optional flag:            shows options

For example, KNeXT can obtain all KGML files for Homo sapiens:

.. code:: text

    $ knext get-kgml hsa

The resulting output folder can be used to parse the files:

.. code:: text
      
    $ knext genes kgml_hsa --graphics

The resulting output folder can be used to convert the TSV files and graphics file:

.. code:: text
      
    $ knext convert kegg_gene_network_hsa hsa --graphics kegg_gene_network_hsa

Graphics text file usage

.. code:: text

    import networkx as nx
    import pandas as pd
    import json
    edges = pd.read_csv('edges.tsv', sep = '\t')
    file_header = open('graphics.txt').read()
    pos = json.loads(file_header)
    graph = nx.from_pandas_edgelist(edges, source = 'entry1', target = 'entry2')
    nx.draw(graph, pos = pos)

Inputs
------

KNeXT only accepts KGML files downloaded from `KEGG <https://www.genome.jp/kegg/>`__

The output of which can be used in successive commands.
All input formats *must be* in TSV format.
Column names are mandatory and should not be changed.

Data Frames
'''''''''''

.. csv-table:: Example TSV file with KEGG ID's
	:header: entry1, entry2, type, value, name

	hsa:100271927-98, hsa:22800-12, PPrel, -->, activation
	hsa:100271927-98, hsa:22808-12, PPrel, -->, activation
	hsa:100271927-98, hsa:3265-12, PPrel, -->, activation

.. csv-table:: Example TSV file for uniprot conversion with `--unique` output 
	:escape: `
        :header: entry1, entry2, type, value, name

	Q9Y243-23, O15111-59, PPrel, -->, activation
	Q9Y243-23, Q6GYQ0-240, PPrel`,`PPrel, --``|```,`+p, inhibition`,`phosphorylation
	Q9Y243-23, O14920-59, PPrel, -->, activation

Installation
------------

The current release is :code:`v1.2.1`
Installation is via pip:

.. code:: bash

    $ pip install knext

Repo can be downloaded and installed through poetry__:

.. code:: bash

    $ git clone https://github.com/everest/knext.git
    $ cd knext
    $ poetry shell
    $ poetry install
    $ poetry run knext [get-kgml, genes, mixed, or convert]

.. __: https://python-poetry.org/

Requirements
------------

Requirements are (also see ``pyproject.toml``):

- Python >= 3.9
- typer__
- click__
- requests__
- pandas__
- networkx__
- pytest__
- pathlib__
- pytest__

.. __: https://typer.tiangolo.com/
.. __: https://click.palletsprojects.com/en/8.1.x/
.. __: https://requests.readthedocs.io/en/latest/
.. __: https://pandas.pydata.org/
.. __: https://networkx.org/
.. __: https://docs.pytest.org/en/7.2.x/
.. __: https://pathlib.readthedocs.io/en/pep428/
.. __: https://docs.pytest.org/en/7.1.x/contents.html
