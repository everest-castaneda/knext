
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining mixed pathway output
"""

import re
import json
import typer
import numpy as np
import pandas as pd
from pathlib import Path
import urllib.request as request
from itertools import combinations
import xml.etree.ElementTree as ET
from collections import defaultdict

def conv_dict(root):
    # This dictionary is the default option
    # Only compounds are unique
    entry_dict = defaultdict(list)
    for entries in root.findall('entry'):
        for key, items in entries.attrib.items():
            entry_dict[key].append(items)  
       
    entry_id=[]
    entry_name=[]
    for key, items in entry_dict.items():
        if key == 'id':
            for i in items:
                entry_id.append(i)
        if key == 'name':
            for i in items:
                entry_name.append(i)

    unique_compound = []
    for i in range(0, len(entry_id)):
        c = [name + '-' + entry_id[i] if name.startswith('cpd:') or name == 'undefined' else name for name in entry_name[i].split()]
        unique_compound.append(' '.join(c))

    conversion_dictionary = dict(zip(entry_id, unique_compound))
    return conversion_dictionary

def conv_dict_unique(root):
    # This dictionary is the unique version
    # Every item is unique to reveal subgraphs
    entry_dict = defaultdict(list)
    for entries in root.findall('entry'):
        for key, items in entries.attrib.items():
            entry_dict[key].append(items)  
       
    entry_id=[]
    entry_name=[]
    for key, items in entry_dict.items():
        if key == 'id':
            for i in items:
                entry_id.append(i)
        if key == 'name':
            for i in items:
                entry_name.append(i)
   
    unique_name = []
    for i in range(0, len(entry_id)):
        e = [name + '-' + entry_id[i] for name in entry_name[i].split()]
        unique_name.append(' '.join(e))

    conversion_dictionary = dict(zip(entry_id, unique_name))
    return conversion_dictionary


def graphics_dict(root):
    entry_dict = defaultdict(list)
    for entry in root.findall('entry'):
        for graphics in entry.find('graphics').items():
            if graphics[0] == 'x' or graphics[0] == 'y':
                entry_dict[entry.attrib['id']].append(int(graphics[1]))
    graphics_dict = {}
    for key, items in entry_dict.items():
        graphics_dict[key] = tuple(items)
    return graphics_dict

def names_dict(root, organism, conversion_dictionary):
    # d = conv_dict_unique(root)
    # d = conv_dict(root)
    e1 = []
    e2 = []
    for entry in root.findall('relation'):
        e1.append(entry.get('entry1'))
        e2.append(entry.get('entry2'))
    e = e1 + e2
    e_list = [conversion_dictionary[entry].split(' ') for entry in e]
    e_list1 = [l for sublist in e_list for l in sublist]
    e_conv = set(e_list1)
    dd = {}
    for n in e_conv:
        # Uses organism code since there are pathways, undefined, and others that will cause
        # an error if used here
        if n.startswith(organism):
            # Remove, if necessary, any terminal modifiers to avoid an error in api call
            n4url = re.sub(r'-[0-9]+', '', n)
            # Uses find to get gene info since other api tools give error
            url = 'https://rest.kegg.jp/find/genes/%s/'
            response = request.urlopen(url % n4url).read().decode('utf-8')
            split_response = response.split('\n')
            # Find only the query gene if given back several accessions that are similar
            s = filter(lambda x: x.startswith(n4url + '\t'), split_response)
            # Adds to dictionary the end entry, which is the written out name
            dd[n] = re.sub('^ ', '', list(s)[0].split(';')[1])
        # Only obtains compounds
        elif n.startswith('cpd:'):
            # Remove terminal modifiers, which are always added to compounds
            # unless a non-unique mixed pathway is chosen
            n4url = re.sub(r'-[0-9]+', '', n)
            # Uses find to get gene info since other api tools give error
            url = 'https://rest.kegg.jp/find/compound/%s'
            response = request.urlopen(url % n4url).read().decode('utf-8')
            subbed_response = re.sub(r'%s\t' % n4url, '', response)
            # Find only the query gene if given back several accessions that are similar
            split_response = re.sub('^ ', '', subbed_response.strip('\n').split(';')[1])
            # Adds to dictionary the end entry, which is the written out name
            dd[n] = split_response
        elif n.startswith('path:'):
            n4url1 = re.sub(r'-[0-9]+', '', n)
            n4url2 = re.sub(r'path:{}'.format(organism), '', n4url1)
            url = 'https://rest.kegg.jp/find/pathway/%s'
            response = request.urlopen(url % n4url2).read().decode('utf-8').strip('\n').split('\t')
            dd[n] = response[1]
        else:
            dd[n] = np.nan
    return dd

def mixed_file(input_data: str, wd: Path, unique: bool = False, graphics: bool = False, names : bool = False):
    """
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of mixed genes, compounds, and pathways that can be used in graph 
    analysis. If -u/--unique flag is used nodes are returned with terminal 
    modifiers to enhance network visualization or analysis. If the 
    -g/--graphics flag is used, x-y coordinates of pathways are returned, 
    which may be used as positions in NetworkX\'s graph drawing commands.
    """
    file = Path(input_data)
    tree = ET.parse(file)
    root = tree.getroot()
   
    title = root.get('title')
    pathway = root.get('name').replace('path:', '')
    pathway_link = root.get('link')
    if graphics:
        # Graphics
        graphics = graphics_dict(root)
        
        typer.echo(f'Now parsing: {title}...')
           
        e1 = []
        e2 = []
        n = []
        v = []
        for relation in root.findall('relation'):
            e1.append(relation.get('entry1'))
            e2.append(relation.get('entry2'))
            for subtype in relation.findall('subtype'):
                n.append(subtype.get('name'))
                v.append(subtype.get('value'))
           
        d = []
        for relation in root.findall('relation'):
            for subtype in relation:
                d1=relation.attrib
                d2=subtype.attrib
                d3=json.dumps(d1),json.dumps(d2)
                d.append(d3)    

        edgelist=[]
        # Takes json dictionary and removes formatting in order to obtain an edgelist
        # This odd line of code is necessary  because of the various nested and non-nested subtypes
        # This allows for edges to inherit weights from a nested subtype
        for line in d:
            x=line[0].replace("{","").replace("}","").replace('"entry1":',"")\
                .replace('"entry2":',"").replace('"type":',"").replace('"name":',"")\
                .replace('"value":',"").replace('"','').replace(',',"\t").replace(' ', '')
            y=line[1].replace("{","").replace("}","").replace('"name":',"")\
                .replace('"',"").replace('value:',"").replace(",",'\t').replace(' ', '')
            edgelist.append(x+"\t"+y)

        df=pd.DataFrame(edgelist)
       
        if df.empty:
            typer.echo(f'File "{file}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')
        else:
            df=df[0].str.split("\t", expand=True).rename({0: 'entry1',1: 'entry2',
                                                      2: 'types', 3:'name',
                                                      4: 'value'}, axis='columns')

            entry_dict = defaultdict(list)
            for entries in root.findall('entry'):
                for key, items in entries.attrib.items():
                    entry_dict[key].append(items) 
            
            # Graphics
            df['pos1'] = df['entry1'].map(graphics)
            df['pos2'] = df['entry2'].map(graphics)
           
            entry_id=[]
            entry_name=[]
            entry_type=[]
            for key, items in entry_dict.items():
                if key == 'id':
                    for i in items:
                        entry_id.append(i)
                if key == 'name':
                    for i in items:
                        entry_name.append(i)
                if key == 'type':
                    for i in items:
                        entry_type.append(i)
            
            if unique:
                conversion_dictionary = conv_dict_unique(root)
            else:
                conversion_dictionary = conv_dict(root)
           
            df['entry1'] = df['entry1'].map(conversion_dictionary)
            df['entry2'] = df['entry2'].map(conversion_dictionary)
            df['entry1'] = df['entry1'].astype(str).str.split(' ', expand = False)
            df['entry2'] = df['entry2'].astype(str).str.split(' ', expand = False)
           
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
                network = [tup + (rows['types'], rows['value'], rows['name'], rows['pos1'], rows['pos2']) for tup in cliques]
                edges.append(network)
            # Standard line for obtaining a dataframe of edges
            edges2df = [e for edge in edges for e in edge]
            # Create pandas DF from edges
            df_out = pd.DataFrame.from_records(edges2df, columns = ['entry1', 'entry2', 'type', 'value', 'name', 'pos1', 'pos2'])
            
            # graphics
            pos_dict1 = {}
            pos_dict2 = {}
            for index, rows in df_out.iterrows():
                pos_dict1[rows['entry1']] = rows['pos1']
                pos_dict2[rows['entry2']] = rows['pos2']
            pos = pos_dict1 | pos_dict2
            json_dict = json.dumps(pos)
            with open(wd / '{}_graphics.txt'.format(pathway), 'w') as outfile:
                outfile.write(json_dict)
            
            dft = df_out.groupby(['entry1', 'entry2'])['type'].apply(list).reset_index()
            dfv = df_out.groupby(['entry1', 'entry2'])['value'].apply(list).reset_index()
            dfn = df_out.groupby(['entry1', 'entry2'])['name'].apply(list).reset_index()
            dfx = dft
            dfx['value'] = dfv['value'].agg(','.join)
            dfx['name'] = dfn['name'].agg(','.join)
            dfx['type'] = dft['type'].agg(','.join)
            # Ensures independently parsed cliques overwrite the cliques, which inherited neighbor weights
            xdf = pd.concat([dfx, cliquedf]).drop_duplicates(subset = ['entry1', 'entry2'], keep = 'last')
            # Removes unneccessary extra "OR" edges connecting to each other from the final dataframe
            # Comment out and remove the 1 from the rest of the dataframes if you want to see their 
            # interaction in the final dataframe, but these are not meant to interact
            xdf_out = xdf[xdf.name != 'clique']
            if names:
                names_dictionary = names_dict(root, root.get('org'), conversion_dictionary)
                xdf_out['entry1_name'] = xdf_out.entry1.map(names_dictionary)
                xdf_out['entry2_name'] = xdf_out.entry2.map(names_dictionary)
                # Cleans up the dataframe for the entry name to be closer
                # to the entry accession code
                xdf_out.insert(1, 'entry1_name', xdf_out.pop('entry1_name'))
                xdf_out.insert(3, 'entry2_name', xdf_out.pop('entry2_name'))
                xdf_out.to_csv(wd / '{}.tsv'.format(pathway), sep = '\t', index = False)
            else:
                xdf_out.to_csv(wd / '{}.tsv'.format(pathway), sep = '\t', index = False)
    else:
        typer.echo(f'Now parsing: {title}...')
           
        e1 = []
        e2 = []
        n = []
        v = []
        for relation in root.findall('relation'):
            e1.append(relation.get('entry1'))
            e2.append(relation.get('entry2'))
            for subtype in relation.findall('subtype'):
                n.append(subtype.get('name'))
                v.append(subtype.get('value'))
           
        d = []
        for relation in root.findall('relation'):
            for subtype in relation:
                d1=relation.attrib
                d2=subtype.attrib
                d3=json.dumps(d1),json.dumps(d2)
                d.append(d3)    

        edgelist=[]
        # Takes json dictionary and removes formatting in order to obtain an edgelist
        # This odd line of code is necessary  because of the various nested and non-nested subtypes
        # This allows for edges to inherit weights from a nested subtype
        for line in d:
            x=line[0].replace("{","").replace("}","").replace('"entry1":',"")\
                .replace('"entry2":',"").replace('"type":',"").replace('"name":',"")\
                .replace('"value":',"").replace('"','').replace(',',"\t").replace(' ', '')
            y=line[1].replace("{","").replace("}","").replace('"name":',"")\
                .replace('"',"").replace('value:',"").replace(",",'\t').replace(' ', '')
            edgelist.append(x+"\t"+y)

        df=pd.DataFrame(edgelist)
       
        if df.empty:
            typer.echo(f'File "{file}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')
        else:
            df=df[0].str.split("\t", expand=True).rename({0: 'entry1',1: 'entry2',
                                                      2: 'types', 3:'name',
                                                      4: 'value'}, axis='columns')

            entry_dict = defaultdict(list)
            for entries in root.findall('entry'):
                for key, items in entries.attrib.items():
                    entry_dict[key].append(items)
           
            entry_id=[]
            entry_name=[]
            entry_type=[]
            for key, items in entry_dict.items():
                if key == 'id':
                    for i in items:
                        entry_id.append(i)
                if key == 'name':
                    for i in items:
                        entry_name.append(i)
                if key == 'type':
                    for i in items:
                        entry_type.append(i)
            
            if unique:
                conversion_dictionary = conv_dict_unique(root)
            else:
                conversion_dictionary = conv_dict(root)
           
            df['entry1'] = df['entry1'].map(conversion_dictionary)
            df['entry2'] = df['entry2'].map(conversion_dictionary)
            df['entry1'] = df['entry1'].astype(str).str.split(' ', expand = False)
            df['entry2'] = df['entry2'].astype(str).str.split(' ', expand = False)
           
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
                network = [tup + (rows['types'], rows['value'], rows['name']) for tup in cliques]
                edges.append(network)
            # Standard line for obtaining a dataframe of edges
            edges2df = [e for edge in edges for e in edge]
            # Create pandas DF from edges
            df_out = pd.DataFrame.from_records(edges2df, columns = ['entry1', 'entry2', 'type', 'value', 'name'])
            
            dft = df_out.groupby(['entry1', 'entry2'])['type'].apply(list).reset_index()
            dfv = df_out.groupby(['entry1', 'entry2'])['value'].apply(list).reset_index()
            dfn = df_out.groupby(['entry1', 'entry2'])['name'].apply(list).reset_index()
            dfx = dft
            dfx['value'] = dfv['value'].agg(','.join)
            dfx['name'] = dfn['name'].agg(','.join)
            dfx['type'] = dft['type'].agg(','.join)
            # Ensures independently parsed cliques overwrite the cliques, which inherited neighbor weights
            xdf = pd.concat([dfx, cliquedf]).drop_duplicates(subset = ['entry1', 'entry2'], keep = 'last')
            # Removes unneccessary extra "OR" edges connecting to each other from the final dataframe
            # Comment out and remove the 1 from the rest of the dataframes if you want to see their 
            # interaction in the final dataframe, but these are not meant to interact
            xdf_out = xdf[xdf.name != 'clique']
            if names:
                names_dictionary = names_dict(root, root.get('org'), conversion_dictionary)
                xdf_out['entry1_name'] = xdf_out.entry1.map(names_dictionary)
                xdf_out['entry2_name'] = xdf_out.entry2.map(names_dictionary)
                # Cleans up the dataframe for the entry name to be closer
                # to the entry accession code
                xdf_out.insert(1, 'entry1_name', xdf_out.pop('entry1_name'))
                xdf_out.insert(3, 'entry2_name', xdf_out.pop('entry2_name'))
                xdf_out.to_csv(wd / '{}.tsv'.format(pathway), sep = '\t', index = False)
            else:
                xdf_out.to_csv(wd / '{}.tsv'.format(pathway), sep = '\t', index = False)

def mixed_folder(input_data: str, wd: Path, unique: bool = False, graphics: bool = False, names: bool = False):
    for file in Path(input_data).glob('*.xml'):
        tree = ET.parse(file)
        root = tree.getroot()
   
        title = root.get('title')
        pathway = root.get('name').replace('path:', '')
        species = root.get('org')
        pathway_link = root.get('link')
        if graphics:
            # Graphics
            graphics = graphics_dict(root)
            
            output_path = wd / 'kegg_mixed_network_{}'.format(species)
            output_path.mkdir(exist_ok = True)
            
            typer.echo(f'Now parsing: {title}...')
           
            e1 = []
            e2 = []
            n = []
            v = []
            for relation in root.findall('relation'):
                e1.append(relation.get('entry1'))
                e2.append(relation.get('entry2'))
                for subtype in relation.findall('subtype'):
                    n.append(subtype.get('name'))
                    v.append(subtype.get('value'))
           
            d = []
            for relation in root.findall('relation'):
                for subtype in relation:
                    d1=relation.attrib
                    d2=subtype.attrib
                    d3=json.dumps(d1),json.dumps(d2)
                    d.append(d3)    

            edgelist=[]
            # Takes json dictionary and removes formatting in order to obtain an edgelist
            # This odd line of code is necessary  because of the various nested and non-nested subtypes
            # This allows for edges to inherit weights from a nested subtype
            for line in d:
                x=line[0].replace("{","").replace("}","").replace('"entry1":',"")\
                .replace('"entry2":',"").replace('"type":',"").replace('"name":',"")\
                .replace('"value":',"").replace('"','').replace(',',"\t").replace(' ', '')
                y=line[1].replace("{","").replace("}","").replace('"name":',"")\
                .replace('"',"").replace('value:',"").replace(",",'\t').replace(' ', '')
                edgelist.append(x+"\t"+y)

            df=pd.DataFrame(edgelist)
       
            if df.empty:
                typer.echo(f'File "{file}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')
            else:
                df=df[0].str.split("\t", expand=True).rename({0: 'entry1',1: 'entry2',
                                                      2: 'types', 3:'name',
                                                      4: 'value'}, axis='columns')

                entry_dict = defaultdict(list)
                for entries in root.findall('entry'):
                    for key, items in entries.attrib.items():
                        entry_dict[key].append(items)  
                
                # Graphics
                df['pos1'] = df['entry1'].map(graphics)
                df['pos2'] = df['entry2'].map(graphics)
           
                entry_id=[]
                entry_name=[]
                entry_type=[]
                for key, items in entry_dict.items():
                    if key == 'id':
                        for i in items:
                            entry_id.append(i)
                    if key == 'name':
                        for i in items:
                            entry_name.append(i)
                    if key == 'type':
                        for i in items:
                            entry_type.append(i)
                
                if unique:
                    conversion_dictionary = conv_dict_unique(root)
                else:
                    conversion_dictionary = conv_dict(root)
           
                df['entry1'] = df['entry1'].map(conversion_dictionary)
                df['entry2'] = df['entry2'].map(conversion_dictionary)
                df['entry1'] = df['entry1'].astype(str).str.split(' ', expand = False)
                df['entry2'] = df['entry2'].astype(str).str.split(' ', expand = False)
           
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
                    network = [tup + (rows['types'], rows['value'], rows['name'], rows['pos1'], rows['pos2']) for tup in cliques]
                    edges.append(network)
            
                # This removes edges which contain +p (phosphorylation) these oftentimes
                # overwrite important weight attributes while providing no vital information
                # edges2df = [e for edge in edges for e in edge if '+p' not in e and '-p' not in e]
            
                # Standard line for obtaining a dataframe of edges
                edges2df = [e for edge in edges for e in edge]
                # Create pandas DF from edges
                df_out = pd.DataFrame.from_records(edges2df, columns = ['entry1', 'entry2', 'type', 'value', 'name', 'pos1', 'pos2'])
                
                # Graphics
                pos_dict1 = {}
                pos_dict2 = {}
                for index, rows in df_out.iterrows():
                    pos_dict1[rows['entry1']] = rows['pos1']
                    pos_dict2[rows['entry2']] = rows['pos2']
                pos = pos_dict1 | pos_dict2
                json_dict = json.dumps(pos)
                with open(output_path / '{}_graphics.txt'.format(pathway), 'w') as outfile:
                    outfile.write(json_dict)
                
                
                dft = df_out.groupby(['entry1', 'entry2'])['type'].apply(list).reset_index()
                dfv = df_out.groupby(['entry1', 'entry2'])['value'].apply(list).reset_index()
                dfn = df_out.groupby(['entry1', 'entry2'])['name'].apply(list).reset_index()
                dfx = dft
                dfx['value'] = dfv['value'].agg(','.join)
                dfx['name'] = dfn['name'].agg(','.join)
                dfx['type'] = dft['type'].agg(','.join)
                # Ensures independently parsed cliques overwrite the cliques, which inherited neighbor weights
                xdf = pd.concat([dfx, cliquedf]).drop_duplicates(subset = ['entry1', 'entry2'], keep = 'last')
                # Removes unneccessary extra "OR" edges connecting to each other from the final dataframe
                # Comment out and remove the 1 from the rest of the dataframes if you want to see their 
                # interaction in the final dataframe, but these are not meant to interact
                xdf_out = xdf[xdf.name != 'clique']
                if names:
                    names_dictionary = names_dict(root, root.get('org'), conversion_dictionary)
                    xdf_out['entry1_name'] = xdf_out.entry1.map(names_dictionary)
                    xdf_out['entry2_name'] = xdf_out.entry2.map(names_dictionary)
                    # Cleans up the dataframe for the entry name to be closer
                    # to the entry accession code
                    xdf_out.insert(1, 'entry1_name', xdf_out.pop('entry1_name'))
                    xdf_out.insert(3, 'entry2_name', xdf_out.pop('entry2_name'))
                    xdf_out.to_csv(output_path / '{}.tsv'.format(pathway), sep = '\t', index = False)
                else:
                    xdf_out.to_csv(output_path / '{}.tsv'.format(pathway), sep = '\t', index = False)
        else:
            output_path = wd / 'kegg_mixed_network_{}'.format(species)
            output_path.mkdir(exist_ok = True)
            
            typer.echo(f'Now parsing: {title}...')
           
            e1 = []
            e2 = []
            n = []
            v = []
            for relation in root.findall('relation'):
                e1.append(relation.get('entry1'))
                e2.append(relation.get('entry2'))
                for subtype in relation.findall('subtype'):
                    n.append(subtype.get('name'))
                    v.append(subtype.get('value'))
           
            d = []
            for relation in root.findall('relation'):
                for subtype in relation:
                    d1=relation.attrib
                    d2=subtype.attrib
                    d3=json.dumps(d1),json.dumps(d2)
                    d.append(d3)    

            edgelist=[]
            # Takes json dictionary and removes formatting in order to obtain an edgelist
            # This odd line of code is necessary  because of the various nested and non-nested subtypes
            # This allows for edges to inherit weights from a nested subtype
            for line in d:
                x=line[0].replace("{","").replace("}","").replace('"entry1":',"")\
                .replace('"entry2":',"").replace('"type":',"").replace('"name":',"")\
                .replace('"value":',"").replace('"','').replace(',',"\t").replace(' ', '')
                y=line[1].replace("{","").replace("}","").replace('"name":',"")\
                .replace('"',"").replace('value:',"").replace(",",'\t').replace(' ', '')
                edgelist.append(x+"\t"+y)

            df=pd.DataFrame(edgelist)
       
            if df.empty:
                typer.echo(f'File "{file}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')
            else:
                df=df[0].str.split("\t", expand=True).rename({0: 'entry1',1: 'entry2',
                                                      2: 'types', 3:'name',
                                                      4: 'value'}, axis='columns')

                entry_dict = defaultdict(list)
                for entries in root.findall('entry'):
                    for key, items in entries.attrib.items():
                        entry_dict[key].append(items)
           
                entry_id=[]
                entry_name=[]
                entry_type=[]
                for key, items in entry_dict.items():
                    if key == 'id':
                        for i in items:
                            entry_id.append(i)
                    if key == 'name':
                        for i in items:
                            entry_name.append(i)
                    if key == 'type':
                        for i in items:
                            entry_type.append(i)
                
                if unique:
                    conversion_dictionary = conv_dict_unique(root)
                else:
                    conversion_dictionary = conv_dict(root)
           
                df['entry1'] = df['entry1'].map(conversion_dictionary)
                df['entry2'] = df['entry2'].map(conversion_dictionary)
                df['entry1'] = df['entry1'].astype(str).str.split(' ', expand = False)
                df['entry2'] = df['entry2'].astype(str).str.split(' ', expand = False)
           
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
                    network = [tup + (rows['types'], rows['value'], rows['name']) for tup in cliques]
                    edges.append(network)
            
                # This removes edges which contain +p (phosphorylation) these oftentimes
                # overwrite important weight attributes while providing no vital information
                # edges2df = [e for edge in edges for e in edge if '+p' not in e and '-p' not in e]
            
                # Standard line for obtaining a dataframe of edges
                edges2df = [e for edge in edges for e in edge]
                # Create pandas DF from edges
                df_out = pd.DataFrame.from_records(edges2df, columns = ['entry1', 'entry2', 'type', 'value', 'name'])
                dft = df_out.groupby(['entry1', 'entry2'])['type'].apply(list).reset_index()
                dfv = df_out.groupby(['entry1', 'entry2'])['value'].apply(list).reset_index()
                dfn = df_out.groupby(['entry1', 'entry2'])['name'].apply(list).reset_index()
                dfx = dft
                dfx['value'] = dfv['value'].agg(','.join)
                dfx['name'] = dfn['name'].agg(','.join)
                dfx['type'] = dft['type'].agg(','.join)
                # Ensures independently parsed cliques overwrite the cliques, which inherited neighbor weights
                xdf = pd.concat([dfx, cliquedf]).drop_duplicates(subset = ['entry1', 'entry2'], keep = 'last')
                # Removes unneccessary extra "OR" edges connecting to each other from the final dataframe
                # Comment out and remove the 1 from the rest of the dataframes if you want to see their 
                # interaction in the final dataframe, but these are not meant to interact
                xdf_out = xdf[xdf.name != 'clique']
                if names:
                    names_dictionary = names_dict(root, root.get('org'), conversion_dictionary)
                    xdf_out['entry1_name'] = xdf_out.entry1.map(names_dictionary)
                    xdf_out['entry2_name'] = xdf_out.entry2.map(names_dictionary)
                    # Cleans up the dataframe for the entry name to be closer
                    # to the entry accession code
                    xdf_out.insert(1, 'entry1_name', xdf_out.pop('entry1_name'))
                    xdf_out.insert(3, 'entry2_name', xdf_out.pop('entry2_name'))
                    xdf_out.to_csv(output_path / '{}.tsv'.format(pathway), sep = '\t', index = False)
                else:
                    xdf_out.to_csv(output_path / '{}.tsv'.format(pathway), sep = '\t', index = False)







