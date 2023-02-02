#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for converting pathway TSV files into UniProt and NCBI IDs
"""

import re
import json
import typer
import urllib.request as request
import pandas as pd
from pathlib import Path
from typing import Optional

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

@app.command()
def file(file: typer.FileText, species, uniprot: bool = False, unique: bool = False, graphics: Optional[Path] = typer.Option(None, exists = True, file_okay = True, dir_okay = False)):
    if uniprot and unique:
        typer.echo(f'Now converting {file.name} to UniProt IDs...')
        df = pd.read_csv(file, delimiter = '\t')
        conv = UP(species)
        df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
        df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
        df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
        df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
        df['entry1'] = df['entry1'].map(conv)
        df['entry2'] = df['entry2'].map(conv)
        df1 = df.dropna()
        for index, row in df1.iterrows():
            df1.loc[index, 'entry1'] = [x.replace('up:', '') + row['match1'] for x in row['entry1']]
            df1.loc[index, 'entry2'] = [y.replace('up:', '') + row['match2'] for y in row['entry2']]
        df_out = df1.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
        df_out.drop(['match1', 'match2'], axis = 1, inplace = True)
        df_out.to_csv('up_%s' % file.name, sep = '\t', index = False)
        if graphics:
            typer.echo(f'Graphics file given! Now converting {graphics.name} to UniProt IDs...')
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
            with open('up_{}'.format(graphics.name), 'w') as outfile:
                outfile.write(json.dumps(up_dict))
    elif not uniprot and unique:
        typer.echo(f'Now converting {file.name} to NCBI IDs...')
        df = pd.read_csv(file, delimiter = '\t')
        conversion = NCBI(species)
        df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
        df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
        df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
        df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
        df['entry1'] = df['entry1'].map(conversion)
        df['entry2'] = df['entry2'].map(conversion)
        df1 = df.dropna()
        df1['entry1'] = df1['entry1'].str.replace('ncbi-geneid:', '') + df1['match1']
        df1['entry2'] = df1['entry2'].str.replace('ncbi-geneid:', '') + df1['match2']
        df2 = df1.drop(['match1', 'match2'], axis = 1) 
        df2.to_csv('ncbi_%s' % file.name, sep = '\t', index = False)
        if graphics:
            typer.echo(f'Graphics file given. Now converting {graphics.name}...')
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
            with open('ncbi_{}'.format(graphics.name), 'w') as outfile:
                outfile.write(json.dumps(ncbi_dict))
    elif uniprot and not unique:
        typer.echo(f'Now converting {file.name} to UniProt IDs...')
        df = pd.read_csv(file, delimiter = '\t')
        conv = UP(species)
        df['entry1'] = df['entry1'].map(conv)
        df['entry2'] = df['entry2'].map(conv)
        df=df.dropna()
        df_out = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
        df_out['entry1'] = df_out['entry1'].str.replace('up:', '')
        df_out['entry2'] = df_out['entry2'].str.replace('up:', '')
        df_out.to_csv('up_%s' % file.name, sep = '\t', index = False)
        if graphics:
            typer.echo(f'Graphics folder given! Now converting {graphics.name} to UniProt IDs...')
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
            with open('up_%s' % graphics.name, 'w') as outfile:
                outfile.write(json.dumps(up_dict))
    elif not uniprot and not unique:
        typer.echo(f'Now converting {file.name} to NCBI IDs...')
        df = pd.read_csv(file, delimiter = '\t')
        conversion = NCBI(species)
        df['entry1'] = df['entry1'].map(conversion).str.replace('ncbi-geneid:', '')
        df['entry2'] = df['entry2'].map(conversion).str.replace('ncbi-geneid:', '')
        df1 = df.dropna()
        df1.to_csv('ncbi_%s' % file.name, sep = '\t', index = False)
        if graphics:
            typer.echo(f'Graphics folder given! Now converting {graphics} to NCBI IDs...')
            pos = open(graphics)
            d = json.loads(pos.read())
            ncbi_dict = {}
            for key, items in d.items():
                try:
                    ncbi_dict[conversion[key].replace('ncbi-geneid:', '')] = items
                except KeyError:
                    pass
            with open('ncbi_%s' % graphics.name, 'w') as outfile:
                outfile.write(json.dumps(ncbi_dict))
        
@app.command()
def folder(fp: Path, species: str, uniprot: bool = False, unique: bool = False, graphics: Optional[Path] = typer.Option(None, exists = True, file_okay = False, dir_okay = True)):
    wd = Path.cwd()
    if uniprot and not unique:
        output_path = wd / 'uniprot_{}'.format(fp)
        output_path.mkdir(exist_ok = True)
        for file in Path(fp).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to UniProt IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conv = UP(species)
            df['entry1'] = df['entry1'].map(conv)
            df['entry2'] = df['entry2'].map(conv)
            df = df.dropna()
            df_out = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
            df_out['entry1'] = df_out['entry1'].str.replace('up:', '')
            df_out['entry2'] = df_out['entry2'].str.replace('up:', '')
            df_out.to_csv(output_path / 'up_{}'.format(file.name), sep = '\t', index = False)
        if graphics and graphics.is_dir():
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
                with open(output_path /'up_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(up_dict))
    elif not uniprot and not unique:
        output_path = wd / 'ncbi_{}'.format(fp)
        output_path.mkdir(exist_ok = True)
        for file in Path(fp).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to NCBI IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conversion = NCBI(species)
            df['entry1'] = df['entry1'].map(conversion).str.replace('ncbi-geneid:', '')
            df['entry2'] = df['entry2'].map(conversion).str.replace('ncbi-geneid:', '')
            df1 = df.dropna()
            df1.to_csv(output_path / 'ncbi_{}'.format(file.name), sep = '\t', index = False)
        if graphics and graphics.is_dir():
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
                with open(output_path / 'ncbi_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(ncbi_dict))
    elif uniprot and unique:
        output_path = wd / 'uniprot_{}'.format(fp)
        output_path.mkdir(exist_ok = True)
        for file in Path(fp).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to UniProt IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conv = UP(species)
            df['match1'] = df['entry1'].str.extract('(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract('(-[0-9]+)')
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry1'] = df['entry1'].map(conv)
            df['entry2'] = df['entry2'].map(conv)
            df1 = df.dropna()
            for index, row in df1.iterrows():
                df1.loc[index, 'entry1'] = [x.replace('up:', '') + row['match1'] for x in row['entry1']]
                df1.loc[index, 'entry2'] = [y.replace('up:', '') + row['match2'] for y in row['entry2']]
            df_out = df1.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
            df_out.drop(['match1', 'match2'], axis = 1, inplace = True)
            df_out.to_csv(output_path / 'up_{}'.format(file.name), sep = '\t', index = False)
        if graphics and graphics.is_dir():
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
                with open(output_path / 'up_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(up_dict))
    elif not uniprot and unique:
        output_path = wd / 'ncbi_{}'.format(fp)
        output_path.mkdir(exist_ok = True)
        for file in Path(fp).glob('*.tsv'):
            typer.echo(f'Now converting {file.name} to NCBI IDs...')
            df = pd.read_csv(file, delimiter = '\t')
            conversion = NCBI(species)
            df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry1'] = df['entry1'].map(conversion)
            df['entry2'] = df['entry2'].map(conversion)
            df1 = df.dropna()
            df1['entry1'] = df1['entry1'].str.replace('ncbi-geneid:', '') + df1['match1']
            df1['entry2'] = df1['entry2'].str.replace('ncbi-geneid:', '') + df1['match2']
            df2 = df1.drop(['match1', 'match2'], axis = 1) 
            df2.to_csv(output_path / 'ncbi_{}'.format(file.name), sep = '\t', index = False)
        if graphics and graphics.is_dir():
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
                with open(output_path / 'ncbi_{}'.format(f.name), 'w') as outfile:
                    outfile.write(json.dumps(ncbi_dict))
                    
if __name__ == '__main__':
    app()
    

