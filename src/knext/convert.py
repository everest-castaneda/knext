#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for converting pathway TSV files into UniProt and NCBI IDs
"""

import re
import json
import typer
import pandas as pd
from pathlib import Path
import urllib.request as request

pd.options.mode.chained_assignment = None
app = typer.Typer()

def UP(species):
    url = 'http://rest.kegg.jp/conv/%s/uniprot'
    response = request.urlopen(url % species).read().decode('utf-8')
    response = response.rstrip().rsplit('\n')
    entrez = []
    uniprot = []
    for resp in response:
        uniprot.append(resp.rsplit()[0])
        entrez.append(resp.rsplit()[1])
    d = {}
    for key, value in zip(entrez, uniprot):
        if key not in d:
            d[key] = [value]
        else:
            d[key].append(value)
    return d

def NCBI(species):
    url = 'http://rest.kegg.jp/conv/%s/ncbi-geneid'
    response = request.urlopen(url % species).read().decode('utf-8')
    response = response.rstrip().rsplit('\n')
    ncbi = []
    kegg = []
    for resp in response:
        ncbi.append(resp.rsplit()[0])
        kegg.append(resp.rsplit()[1])
    d = {}
    for key, value in zip(kegg, ncbi):
        if key not in d:
            d[key] = value
        else:
            d[key].append(value)
    return d

def convert_file(species, input_data, wd: Path, graphics = None, uniprot: bool = False, unique: bool = False):
    """
    Converts a file or folder of parsed genes or mixed pathways from KEGG IDs to 
    NCBI gene IDs or UniProt IDs. Default is NCBI gene IDs unless the
    -up/--uniprot flag is used. The -g/--graphics flag converts the graphics file
    or folder into the accompanying IDs. Do not mix/match! Use the -u/--unique 
    flag only if the input is unique, presence of terminal modifiers. 
    Use the -g/--graphics flag with folders if the input for the conversion is 
    a folder and vice-versa.
    """
    if Path(input_data).is_file():
        file = Path(input_data)
        if uniprot and unique:
            typer.echo(f'Now converting {file.name} to UniProt IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conv = UP(species)
            # Extract the terminal modifiers and create a new column
            # This enables the re-addition of the modifiers at the
            # end of the function.
            df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
            # Remove the terminal modifier so that the IDs map properly
            # to the KEGG API call
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
            # Map to convert KEGG IDs to UniProt IDs. Note lists are returned
            # for some conversions.
            df['entry1_up'] = df['entry1'].map(conv)
            df['entry2_up'] = df['entry2'].map(conv)
            # Fills nans with entries from original columns
            df['entry1'] = df['entry1_up'].fillna(df['entry1']).apply(lambda x: [x] if type(x) is not list else x)
            df['entry2'] = df['entry2_up'].fillna(df['entry2']).apply(lambda x: [x] if type(x) is not list else x)
            # Drop the extra column as it's all now in entry1/2 columns
            df1 = df.drop(['entry1_up', 'entry2_up'], axis = 1)
            # Drop any IDs that did not convert. (old version)
            # df1 = df.dropna()
            # Joins the comma seperated list to a string to avoid iterrating
            df1['entry1'] = [','.join(map(str, l)) for l in df1['entry1']]
            df1['entry2'] = [','.join(map(str, l)) for l in df1['entry2']]
            # Quick replace of all up: modifiers at the beginning of each string (old method)
            # df1 = df1.replace(to_replace = 'up:', value = '', regex = True)
            # Quick seperation back into a list while avoiding cpd:, path:, and undefined 
            # Also removes up: modifier to avoid breaking explode
            df1['entry1'] = df1['entry1'].apply(lambda x: re.findall(r'[a-zA-z0-9]+', x.replace('up:', '')) if x.startswith('up:') else [x])
            df1['entry2'] = df1['entry2'].apply(lambda x: re.findall(r'[a-zA-z0-9]+', x.replace('up:', '')) if x.startswith('up:') else [x])
            df_out = df1.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
            df_out['entry1'] = df_out['entry1'] + df_out['match1']
            df_out['entry2'] = df_out['entry2'] + df_out['match2']
            df_out.drop(['match1', 'match2'], axis = 1, inplace = True)
            # Finally, remove all rows with 'hsa:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'hsa:' for the for loop
            df_out = df_out[~df_out['entry1'].astype(str).str.startswith(species)]
            df_out = df_out[~df_out['entry2'].astype(str).str.startswith(species)]
            df_out.to_csv(wd / 'up_{}'.format(file.name), sep = '\t', index = False)
            if graphics != None and Path(graphics).is_file():
                typer.echo(f'Graphics file given! Now converting {Path(graphics).name} to UniProt IDs...')
                pos = open(graphics)
                d = json.loads(pos.read())
                up_dict = {}
                for key, items in d.items():
                    pattern = re.search(r'(-[0-9]+)', key)
                    new_keys = re.sub(r'(-[0-9]+)', '', key)
                    try:
                        uplist = conv[new_keys]
                        new_up = [up + pattern.group() for up in uplist]
                        up_dict[str(new_up)] = items
                        for up in new_up:
                            up_dict[up.replace('up:', '')] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        up_dict.update({key: items})
                with open(wd / 'up_{}'.format(Path(graphics).name), 'w') as outfile:
                    outfile.write(json.dumps(up_dict))
        elif not uniprot and unique:
            typer.echo(f'Now converting {file.name} to NCBI IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conversion = NCBI(species)
            df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry1_ncbi'] = df['entry1'].map(conversion)
            df['entry2_ncbi'] = df['entry2'].map(conversion)
            # Fills nans with entries from original columns
            df['entry1'] = df['entry1_ncbi'].fillna(df['entry1'])
            df['entry2'] = df['entry2_ncbi'].fillna(df['entry2'])
            # Drop the extra column as it's all now in entry1/2 columns
            df = df.drop(['entry1_ncbi', 'entry2_ncbi'], axis = 1)
            # Finally, remove all rows with 'species:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'species:' for the for loop
            df1 = df[~df['entry1'].astype(str).str.startswith(species)]
            df1 = df[~df['entry2'].astype(str).str.startswith(species)]
            # df1 = df.dropna()
            df1['entry1'] = df1['entry1'].str.replace('ncbi-geneid:', '') + df1['match1']
            df1['entry2'] = df1['entry2'].str.replace('ncbi-geneid:', '') + df1['match2']
            df_out = df1.drop(['match1', 'match2'], axis = 1) 
            df_out.to_csv(wd / 'ncbi_{}'.format(file.name), sep = '\t', index = False)
            if graphics != None and Path(graphics).is_file():
                typer.echo(f'Graphics file given. Now converting {Path(graphics).name}...')
                pos = open(graphics)
                d = json.loads(pos.read())
                ncbi_dict = {}
                for key, items in d.items():
                    pattern = re.search(r'(-[0-9]+)', key)
                    new_keys = re.sub(r'(-[0-9]+)', '', key)
                    try:
                        ncbi_dict[conversion[new_keys].replace('ncbi-geneid:', '') + pattern.group()] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        ncbi_dict.update({key: items})
                with open(wd / 'ncbi_{}'.format(Path(graphics).name), 'w') as outfile:
                    outfile.write(json.dumps(ncbi_dict))
        elif uniprot and not unique:
            typer.echo(f'Now converting {file.name} to UniProt IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conv = UP(species)
            df['entry1_up'] = df['entry1'].map(conv)
            df['entry2_up'] = df['entry2'].map(conv)
            # Fills nans with entries from original columns
            df['entry1'] = df['entry1_up'].fillna(df['entry1']).apply(lambda x: [x] if type(x) is not list else x)
            df['entry2'] = df['entry2_up'].fillna(df['entry2']).apply(lambda x: [x] if type(x) is not list else x)
            # Drop the extra column as it's all now in entry1/2 columns
            df_out = df.drop(['entry1_up', 'entry2_up'], axis = 1)
            # df=df.dropna()
            df_out = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
            df_out['entry1'] = df_out['entry1'].str.replace('up:', '')
            df_out['entry2'] = df_out['entry2'].str.replace('up:', '')
            # Finally, remove all rows with 'hsa:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'hsa:' for the for loop
            df_out = df_out[~df_out['entry1'].astype(str).str.startswith(species)]
            df_out = df_out[~df_out['entry2'].astype(str).str.startswith(species)]
            df_out.to_csv(wd / 'up_{}'.format(file.name), sep = '\t', index = False)
            if graphics != None and Path(graphics).is_file():
                typer.echo(f'Graphics file given! Now converting {Path(graphics).name} to UniProt IDs...')
                pos = open(graphics)
                d = json.loads(pos.read())
                up_dict = {}
                for key, items in d.items():
                    try:
                        uplist = conv[key]
                        for up in uplist:
                            up_dict[up.replace('up:', '')] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        up_dict.update({key: items})
                with open(wd / 'up_{}'.format(Path(graphics).name), 'w') as outfile:
                    outfile.write(json.dumps(up_dict))
        else:
            typer.echo(f'Now converting {file.name} to NCBI IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conversion = NCBI(species)
            # Creates a new column for conversions, which will likely have nans
            df['entry1_ncbi'] = df['entry1'].map(conversion).str.replace('ncbi-geneid:', '')
            df['entry2_ncbi'] = df['entry2'].map(conversion).str.replace('ncbi-geneid:', '')
            # Fills nans with entries from original columns
            df['entry1'] = df['entry1_ncbi'].fillna(df['entry1'])
            df['entry2'] = df['entry2_ncbi'].fillna(df['entry2'])
            # Drop the extra column as it's all now in entry1/2 columns
            df = df.drop(['entry1_ncbi', 'entry2_ncbi'], axis = 1)
            # Finally, remove all rows with 'species:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'species:' for the for loop
            df1 = df[~df['entry1'].astype(str).str.startswith(species)]
            df1 = df[~df['entry2'].astype(str).str.startswith(species)]
            df1.to_csv(wd / 'ncbi_{}'.format(file.name), sep = '\t', index = False)
            if graphics != None and Path(graphics).is_file():
                typer.echo(f'Graphics file given! Now converting {Path(graphics).name} to NCBI IDs...')
                pos = open(graphics)
                d = json.loads(pos.read())
                ncbi_dict = {}
                for key, items in d.items():
                    try:
                        ncbi_dict[conversion[key].replace('ncbi-geneid:', '')] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        ncbi_dict.update({key: items})
                with open(wd / 'ncbi_{}'.format(Path(graphics).name), 'w') as outfile:
                    outfile.write(json.dumps(ncbi_dict))

def convert_folder(species, input_data, wd: Path, graphics = None, uniprot: bool = False, unique: bool = False):
    if uniprot and not unique:
        output_path = wd / 'up_{}'.format(species)
        output_path.mkdir(exist_ok = True)
        for file in Path(input_data).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to UniProt IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conv = UP(species)
            df['entry1_up'] = df['entry1'].map(conv)
            df['entry2_up'] = df['entry2'].map(conv)
            df['entry1'] = df['entry1_up'].fillna(df['entry1']).apply(lambda x: [x] if type(x) is not list else x)
            df['entry2'] = df['entry2_up'].fillna(df['entry2']).apply(lambda x: [x] if type(x) is not list else x)
            # Drop the extra column as it's all now in entry1/2 columns
            df_out = df.drop(['entry1_up', 'entry2_up'], axis = 1)
            # df = df.dropna()
            df_out = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
            df_out['entry1'] = df_out['entry1'].str.replace('up:', '')
            df_out['entry2'] = df_out['entry2'].str.replace('up:', '')
            # Finally, remove all rows with 'hsa:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'hsa:' for the for loop
            df_out = df_out[~df_out['entry1'].astype(str).str.startswith(species)]
            df_out = df_out[~df_out['entry2'].astype(str).str.startswith(species)]
            df_out.to_csv(output_path / 'up_{}'.format(file.name), sep = '\t', index = False)
        if graphics != None and Path(graphics).is_dir():
            for f in Path(graphics).glob('*.txt'):
                typer.echo(f'Graphics folder given! Now converting {f.name} to UniProt IDs...')
                pos = open(f)
                d = json.loads(pos.read())
                up_dict = {}
                for key, items in d.items():
                    try:
                        uplist = conv[key]
                        for up in uplist:
                            up_dict[up.replace('up:', '')] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        up_dict.update({key: items})
                with open(output_path /'up_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(up_dict))
    elif not uniprot and not unique:
        output_path = wd / 'ncbi_{}'.format(species)
        output_path.mkdir(exist_ok = True)
        for file in Path(input_data).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to NCBI IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conversion = NCBI(species)
            df['entry1_ncbi'] = df['entry1'].map(conversion).str.replace('ncbi-geneid:', '')
            df['entry2_ncbi'] = df['entry2'].map(conversion).str.replace('ncbi-geneid:', '')
            df['entry1'] = df['entry1_ncbi'].fillna(df['entry1'])
            df['entry2'] = df['entry2_ncbi'].fillna(df['entry2'])
            # Drop the extra column as it's all now in entry1/2 columns
            df = df.drop(['entry1_ncbi', 'entry2_ncbi'], axis = 1)
            # Finally, remove all rows with 'species:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'species:' for the for loop
            df1 = df[~df['entry1'].astype(str).str.startswith(species)]
            df1 = df[~df['entry2'].astype(str).str.startswith(species)]
            df1 = df.dropna()
            df1.to_csv(output_path / 'ncbi_{}'.format(file.name), sep = '\t', index = False)
        if graphics != None and Path(graphics).is_dir():
            for f in Path(graphics).glob('*.txt'):
                typer.echo(f'Graphics folder given! Now converting {f.name} to NCBI IDs...')
                pos = open(f)
                d = json.loads(pos.read())
                ncbi_dict = {}
                for key, items in d.items():
                    try:
                        ncbi_dict[conversion[key].replace('ncbi-geneid:', '')] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        ncbi_dict.update({key: items})
                with open(output_path / 'ncbi_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(ncbi_dict))
    elif uniprot and unique:
        output_path = wd / 'up_{}'.format(species)
        output_path.mkdir(exist_ok = True)
        for file in Path(input_data).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to UniProt IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conv = UP(species)
            # Extract the terminal modifiers and create a new column
            # This enables the readdition of the modifiers at the
            # end of the function.
            df['match1'] = df['entry1'].str.extract('(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract('(-[0-9]+)')
            # Remove the terminal modifier so that the IDs map properly
            # to the KEGG API call
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
            # Map to convert KEGG IDs to UniProt IDs. Note lists are returned
            # for some conversions.
            df['entry1_up'] = df['entry1'].map(conv)
            df['entry2_up'] = df['entry2'].map(conv)
            df['entry1'] = df['entry1_up'].fillna(df['entry1']).apply(lambda x: [x] if type(x) is not list else x)
            df['entry2'] = df['entry2_up'].fillna(df['entry2']).apply(lambda x: [x] if type(x) is not list else x)
            # Drop the extra column as it's all now in entry1/2 columns
            df1 = df.drop(['entry1_up', 'entry2_up'], axis = 1)
            # Drop any IDs that did not convert.
            # This avoids errors in the subsequent entry that creates lists
            # df1 = df.dropna()
            # Joins the comma seperated list to a string to avoid iterrating
            df1['entry1'] = [','.join(map(str, l)) for l in df1['entry1']]
            df1['entry2'] = [','.join(map(str, l)) for l in df1['entry2']]
            # Quick replace of all up: modifiers at the beginning of each string (old method)
            # df1 = df1.replace(to_replace = 'up:', value = '', regex = True)
            # Quick seperation back into a list while avoiding cpd:, path:, and undefined 
            # Also removes up: modifier to avoid breaking explode
            df1['entry1'] = df1['entry1'].apply(lambda x: re.findall(r'[a-zA-z0-9]+', x.replace('up:', '')) if x.startswith('up:') else [x])
            df1['entry2'] = df1['entry2'].apply(lambda x: re.findall(r'[a-zA-z0-9]+', x.replace('up:', '')) if x.startswith('up:') else [x])
            # Individualize each entry from a list
            df_out = df1.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
            # Add back the unique modifiers
            df_out['entry1'] = df_out['entry1'] + df_out['match1']
            df_out['entry2'] = df_out['entry2'] + df_out['match2']
            # Drop the unique modifiers columns
            df_out.drop(['match1', 'match2'], axis = 1, inplace = True)
            # Finally, remove all rows with 'hsa:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'hsa:' for the for loop
            df_out = df_out[~df_out['entry1'].astype(str).str.startswith(species)]
            df_out = df_out[~df_out['entry2'].astype(str).str.startswith(species)]
            df_out.to_csv(output_path / 'up_{}'.format(file.name), sep = '\t', index = False)
        if graphics != None and Path(graphics).is_dir():
            for f in Path(graphics).glob('*.txt'):
                typer.echo(f'Graphics folder given! Now converting {f.name}...')
                pos = open(f)
                d = json.loads(pos.read())
                up_dict = {}
                for key, items in d.items():
                    pattern = re.search(r'(-[0-9]+)', key)
                    new_keys = re.sub(r'(-[0-9]+)', '', key)
                    try:
                        uplist = conv[new_keys]
                        new_up = [up + pattern.group() for up in uplist]
                        up_dict[str(new_up)] = items
                        for up in new_up:
                            up_dict[up.replace('up:', '')] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        up_dict.update({key: items})
                with open(output_path / 'up_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(up_dict))
    elif not uniprot and unique:
        output_path = wd / 'ncbi_{}'.format(species)
        output_path.mkdir(exist_ok = True)
        for file in Path(input_data).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to NCBI IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conversion = NCBI(species)
            df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry1_ncbi'] = df['entry1'].map(conversion)
            df['entry2_ncbi'] = df['entry2'].map(conversion)
            df['entry1'] = df['entry1_ncbi'].fillna(df['entry1'])
            df['entry2'] = df['entry2_ncbi'].fillna(df['entry2'])
            # Drop the extra column as it's all now in entry1/2 columns
            df = df.drop(['entry1_ncbi', 'entry2_ncbi'], axis = 1)
            # Finally, remove all rows with 'species:' since this will create misleading files
            # Also clash with the graphics file since it won't include 'species:' for the for loop
            df1 = df[~df['entry1'].astype(str).str.startswith(species)]
            df1 = df[~df['entry2'].astype(str).str.startswith(species)]
            # df1 = df.dropna()
            df1['entry1'] = df1['entry1'].str.replace('ncbi-geneid:', '') + df1['match1']
            df1['entry2'] = df1['entry2'].str.replace('ncbi-geneid:', '') + df1['match2']
            df2 = df1.drop(['match1', 'match2'], axis = 1) 
            df2.to_csv(output_path / 'ncbi_{}'.format(file.name), sep = '\t', index = False)
        if graphics != None and Path(graphics).is_dir():
            for f in Path(graphics).glob('*.txt'):
                typer.echo(f'Graphics folder given! Now converting {f.name} to NCBI IDs...')
                pos = open(f)
                d = json.loads(pos.read())
                ncbi_dict = {}
                for key, items in d.items():
                    pattern = re.search(r'(-[0-9]+)', key)
                    new_keys = re.sub(r'(-[0-9]+)', '', key)
                    try:
                        ncbi_dict[conversion[new_keys].replace('ncbi-geneid:', '') + pattern.group()] = items
                    except KeyError:
                        pass
                for key, items in d.items():
                    if not key.startswith(species):
                        ncbi_dict.update({key: items})
                with open(output_path / 'ncbi_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(ncbi_dict))
                    
                        
                        
            
    
    
    
    