#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining gene-only pathways
"""

import re
import json
import typer
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import urllib.request as request
from itertools import combinations
import xml.etree.ElementTree as ET
from collections import defaultdict
from knext import utils
from knext.utils import FileNotFound




class GenesInteractionParser:

    def __init__(self, input_data: str, wd: Path, mixed:bool = False, unique: bool = False,
                 graphics: bool = False, names: bool = False, verbose: bool = False):
        self.input_data = input_data
        self.wd = wd
        self.mixed = mixed
        self.unique = unique
        self.graphics = graphics
        self.names = names
        self.verbose = verbose

        tree = ET.parse(input_data)
        self.root = tree.getroot()

        self.conversion_dictionary = self._get_conversion_dictionary()
        if self.names:
            self.names_dictionary = self._get_names_dictionary(self.conversion_dictionary)


    def _get_edges(self):
        '''
        This function takes the root of elementTree and returns an edgelist
        '''
        pathway_link = self.root.get('link')

        d = []
        for relation in self.root.findall('relation'):
            for subtype in relation:
                d1=relation.attrib
                d2=subtype.attrib
                d3=json.dumps(d1),json.dumps(d2)
                d.append(d3)

        edgelist=[]
        for line in d:
            x=line[0].replace("{","").replace("}","").replace('"entry1":',"") \
                .replace('"entry2":',"").replace('"type":',"").replace('"name":',"") \
                .replace('"value":',"").replace('"','').replace(',',"\t").replace(' ', '')
            y=line[1].replace("{","").replace("}","").replace('"name":',"") \
                .replace('"',"").replace('value:',"").replace(",",'\t').replace(' ', '')
            edgelist.append(x+"\t"+y)

        df=pd.DataFrame(edgelist)
        if df.empty:
            # throw error if no edges are found
            raise FileNotFound(f'ERROR: File "{self.input_data}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')

        df=df[0].str.split("\t", expand=True).rename({0: 'entry1',1: 'entry2',
                                                      2: 'types', 3:'name',
                                                      4: 'value'}, axis='columns')
        # reorder columns as entry1, entry2, types, value, name
        df = df[['entry1', 'entry2', 'types', 'value', 'name']]
        # convert compound value to kegg id if only relation.type is "compound"
        def apply_conversion(row):
            if row['name'] == 'compound':
                return self.conversion_dictionary.get(row['value'], row['value'])
            else:
                return row['value']
        df['value'] = df.apply(apply_conversion, axis=1)

        # Convert entry1 and entry2 id to kegg id
        df['entry1'] = df['entry1'].map(self.conversion_dictionary)
        df['entry2'] = df['entry2'].map(self.conversion_dictionary)
        # Split the entry1 and entry2 into lists
        # entry1 and entry2 can be a list of genes
        df['entry1'] = df['entry1'].astype(str).str.split(' ', expand = False)
        df['entry2'] = df['entry2'].astype(str).str.split(' ', expand = False)
        return df

    def _get_conversion_dictionary(self):
        if self.unique:
            conversion_dictionary = utils.conv_dict_unique(self.root)
        else:
            conversion_dictionary = utils.conv_dict(self.root)
        return conversion_dictionary

    def _get_names_dictionary(self, conversion_dictionary):
        '''
        Get the names dictionary for the given GenesInteractionParser object
        '''
        names_dictionary = utils.names_dict(self.root, self.root.get('org'), conversion_dictionary)
        return self.names_dictionary


    def _parse_clique(self, df):
        # Parse the cliques seperate so they won't inherit neighbor weights
        # Allows for custom weights so user knows which are customly parsed
        clique_edges = []
        for index, rows in df.iterrows():
            if len(rows['entry1']) > 1:
                cliques = [x for x in combinations(rows['entry1'], 2)]
                cliques1 = [tup + ('type 2', 'undirectional', 'clique') for tup in cliques]
                clique_edges.append(cliques1)
            if len(rows['entry2']) > 1:
                cliques = [x for x in combinations(rows['entry2'], 2)]
                cliques2 = [tup + ('type 2', 'undirectional', 'clique') for tup in cliques]
                clique_edges.append(cliques2)
        clique2df = [e for edge in clique_edges for e in edge]
        cliquedf = pd.DataFrame.from_records(clique2df, columns = ['entry1', 'entry2', 'type', 'value', 'name'])

        edges = []
        for index, rows in df.iterrows():
            lists = rows['entry1'] + rows['entry2']
            # Cliques here inherit neighbor weights, but will be overwritten by above
            cliques = [x for x in combinations(lists, 2)]
            if self.graphics:
                network = [tup + (rows['types'], rows['value'], rows['name'], rows['pos1'], rows['pos2']) for tup in cliques]
            else:
                network = [tup + (rows['types'], rows['value'], rows['name']) for tup in cliques]
            edges.append(network)

        # This removes edges which contain +p (phosphorylation) these oftentimes
        # overwrite important weight attributes while providing no vital information
        # edges2df = [e for edge in edges for e in edge if '+p' not in e and '-p' not in e]

        # Standard line for obtaining a dataframe of edges
        edges2df = [e for edge in edges for e in edge]
        # Create pandas DF from edges
        if self.graphics:
            df_out = pd.DataFrame.from_records(edges2df, columns = ['entry1', 'entry2', 'type', 'value', 'name', 'pos1', 'pos2'])
        else:
            df_out = pd.DataFrame.from_records(edges2df, columns = ['entry1', 'entry2', 'type', 'value', 'name'])
        return cliquedf, df_out

    def _propagate_compounds(self, xdf):
        # These next series of steps are for propagating compounds and undefined nodes
        # Uses networkx
        G = nx.from_pandas_edgelist(xdf, source = 'entry1', target = 'entry2', edge_attr = 'name', create_using = nx.DiGraph())
        new_edges = []
        for node in G.nodes:
            # Gather roots and leaflets
            roots = [n for n, d in G.in_degree() if d == 0]
            leaflets = [n for n, d in G.out_degree() if d == 0]
            # Find compounds or undefined proteins that might need propagation
            if node.startswith('cpd') or node.startswith('undefined'):
                if node in roots or node in leaflets:
                    # Passes root nodes and leaflet nodes
                    # These do not need propagation as they terminate
                    pass
                else:
                    out_edges = [e for e in G.out_edges(node)]
                    in_edges = [e for e in G.in_edges(node)]
                    for i in in_edges:
                        for o in out_edges:
                            if (not i[0].startswith('cpd') and not o[1].startswith('cpd') and
                                    not i[0].startswith('undefined') and
                                    not o[1].startswith('undefined') and
                                    not i[0].startswith('path') and not o[1].startswith('path')):
                                # Simple compound propagation removes compound between two genes
                                # Example: hsa:xxx -> cpd:xxx -> hsa:xxx to hsa:xxx -> hsa:xxx
                                new_edges.append([i[0], o[1], 'CPp', 'Custom', 'compound propagation'])
                            else:
                                # Loops through roots and leaves to find the start and end of each node that is propagated
                                # These pathways have copious amounts of compounds and undefined genes in between the propagated genes
                                # This loop also picks up on non-proteins that are left out if the above simple compound propagation is used
                                for root in roots:
                                    for leaf in leaflets:
                                        if nx.has_path(G, root, node) == True and nx.has_path(G, node, leaf) == True:
                                            # Uses the shortest path between the root/leaf and query node
                                            rpath = [p for p in nx.shortest_path(G, root, node)]
                                            lpath = [p for p in nx.shortest_path(G, node, leaf)]
                                            # Skip paths that have no proteins
                                            if all(node.startswith('cpd') or node.startswith('undefined') or node.startswith('path') for node in rpath) == True or all(node.startswith('cpd') or node.startswith('undefined') or node.startswith('path') for node in lpath) == True:
                                                pass
                                            else:
                                                # Indexes all elements in the shortest path list that is a gene
                                                rindex = [index for index, element in enumerate(rpath) if not element.startswith('cpd') and not element.startswith('undefined') and not element.startswith('path')]
                                                lindex = [index for index, element in enumerate(lpath) if not element.startswith('cpd') and not element.startswith('undefined') and not element.startswith('path')]
                                                # Here we use the maximum of the root path and the minimum of the leaf path to successfully propagate compounds
                                                # Note below example
                                                # [mmu:x, cpd:###, cpd:x] --> max is mmu:x
                                                # [cpd:x, cpd:###, undefined, mmu:y, mmu:###] --> min is mmu:y
                                                # creates edge [mmu:x, mmu:y]
                                                new_edges.append([rpath[max(rindex)], lpath[min(lindex)], 'CPp', 'Custom', 'compound propagation'])
        # Creates a dataframe of the new edges that are a result of compound propagation
        new_edges_df = pd.DataFrame(new_edges, columns = ['entry1', 'entry2', 'type', 'value', 'name'])
        # Concatenates the new edges with the edges from the above (cliques and original parsed edges)
        df0 = pd.concat([xdf, new_edges_df])
        # Drop any duplicated edges
        df1 = df0.drop_duplicates()
        # Removes compounds and undefined as they were propagated and no longer needed
        df2 = df1[(~df1['entry1'].str.startswith('cpd')) & (~df1['entry2'].str.startswith('cpd')) & (~df1['entry1'].str.startswith('undefined')) & (~df1['entry2'].str.startswith('undefined')) & (~df1['entry1'].str.startswith('path')) & (~df1['entry2'].str.startswith('path'))]
        # Removes unneccessary extra "OR" edges connecting to each other from the final dataframe
        # Comment out and remove the 1 from the rest of the dataframes if you want to see their
        # interaction in the final dataframe, but these are not meant to interact
        return df2

    def _replace_with_cliques(self, df, cliquedf, df_out):
        '''
        This function replaces the edges in df with the cliques in cliquedf
        '''
        # aggregating the 'type', 'value', and 'name' columns in df_out for each
        # unique pair of 'entry1' and 'entry2',
        #  joining the aggregated values into a single string,
        # and then merging this with cliquedf while removing any duplicates.
        dft = df_out.groupby(['entry1', 'entry2'])['type'].apply(list).reset_index()
        dfv = df_out.groupby(['entry1', 'entry2'])['value'].apply(list).reset_index()
        dfn = df_out.groupby(['entry1', 'entry2'])['name'].apply(list).reset_index()
        dfx = dft
        dfx['type'] = dft['type'].transform(','.join)
        dfx['value'] = dfv['value'].transform(','.join)
        dfx['name'] = dfn['name'].transform(','.join)
        # Ensures independently parsed cliques overwrite the cliques, which inherited neighbor weights
        xdf = pd.concat([dfx, cliquedf]).drop_duplicates(subset = ['entry1', 'entry2'], keep = 'last')
        return  xdf

    def _add_names(self, df):
        df['entry1_name'] = df.entry1.map(self.names_dictionary)
        df['entry2_name'] = df.entry2.map(self.names_dictionary)
        # Cleans up the dataframe for the entry name to be closer
        # to the entry accession code
        df.insert(1, 'entry1_name', df.pop('entry1_name'))
        df.insert(3, 'entry2_name', df.pop('entry2_name'))
        return df

    def genes_file(self):
        '''
        Converts a folder of KGML files or a single KGML file into a weighted
        edgelist of genes that can be used in graph analysis. If -u/--unique flag
        is used genes are returned with terminal modifiers to enhance network
        visualization or analysis. If the -g/--graphics flag is used, x-y
        coordinates of pathways are returned, which may be used as positions in
        NetworkX\'s graph drawing commands.
        input_data: str: Path to KGML file or folder of KGML files
        wd: Path: Path to working directory
        unique: bool: Flag to return unique gene names
        '''
        title = self.root.get('title')
        pathway = self.root.get('name').replace('path:', '')
        pathway_link = self.root.get('link')

        # Common operations
        if self.verbose:
            typer.echo(typer.style(f'Now parsing: {title}...', fg=typer.colors.GREEN, bold=False))
        df = self._get_edges()

        if self.graphics:
            graphics = utils.graphics_dict(self.root)
            df['pos1'] = df['entry1'].map(graphics)
            df['pos2'] = df['entry2'].map(graphics)

        cliquedf, df_out = self._parse_clique(df)

        if self.graphics:
            self._parse_graphics(df, graphics, df_out)

        xdf = self._replace_with_cliques(df, cliquedf, df_out)




        # Check for compounds or undefined nodes
        has_compounds_or_undefined = not xdf[(xdf['entry1'].str.startswith('cpd:')) | (xdf['entry2'].str.startswith('cpd:')) | (xdf['entry1'].str.startswith('undefined')) | (xdf['entry2'].str.startswith('undefined'))].empty

        if not self.mixed:
            # remove edge with "path" entries
            xdf = xdf[(~xdf['entry1'].str.startswith('path')) &
                      (~xdf['entry2'].str.startswith('path'))]
            if has_compounds_or_undefined:

                xdf = self._propagate_compounds(xdf)
                xdf = xdf[xdf.name != 'clique']
                if self.names:
                    xdf['entry1_name'] = xdf.entry1.map(self.names_dictionary)
        else:
            xdf = xdf[xdf.name != 'clique']
            if self.names:
                xdf = self._add_names(xdf)
        xdf.to_csv(self.wd / '{}.tsv'.format(pathway), sep = '\t', index = False)


    def _parse_graphics(df, graphics, df_out, wd, pathway):
        # Graphics
        pos_dict1 = {}
        pos_dict2 = {}
        for index, rows in df_out.iterrows():
            pos_dict1[rows['entry1']] = rows['pos1']
            pos_dict2[rows['entry2']] = rows['pos2']
        pos = pos_dict1 | pos_dict2
        json_dict = json.dumps(pos)
        with open(wd / '{}_graphics.txt'.format(pathway), 'w') as outfile:
            outfile.write(json_dict)





def genes_parser(input_data: str, wd: Path, mixed:bool = False, unique: bool = False,
                 graphics: bool = False, names: bool = False, verbose: bool = False):
    '''
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of genes that can be used in graph analysis.
    '''
    if Path(input_data).is_dir():
        for file in Path(input_data).glob('*.xml'):
            try:
                gip = GenesInteractionParser(file, wd, mixed=mixed,
                                             unique=unique, graphics=graphics, names=names,
                                             verbose=verbose)
                gip.genes_file()
            except FileNotFound as e:
                typer.echo(typer.style(e.message, fg=typer.colors.RED, bold=True))
                continue
    else:
        gip = GenesInteractionParser(input_data, wd, mixed=mixed,
                                     unique=unique, graphics=graphics, names=names,
                                     verbose=verbose)
        gip.genes_file()