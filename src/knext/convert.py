#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for converting pathway TSV files into UniProt and NCBI IDs
"""

import json
import pathlib
import re
from pathlib import Path

import pandas as pd
import typer

from knext.utils import UP, NCBI, FileNotFound

pd.options.mode.chained_assignment = None
app = typer.Typer()




class Converter:
    def __init__(self, species, input_data, wd: Path, graphics=None,
                 uniprot: bool = False, unique: bool = False, verbose: bool = False):
        self.species = species
        self.input_data = input_data
        self.wd = wd
        self.graphics = graphics
        self.uniprot = uniprot
        self.unique = unique
        if uniprot:
            self.conversion = UP(self.species)
            self.prefix = 'up:'
        else:
            self.conversion = conversion = NCBI(self.species)
            self.prefix = 'ncbi-geneid:'

    def _process_graphics(self):
        # extract the filename part of self.input_data
        graphics_file =pathlib.PurePath(self.graphics, Path(self.input_data).stem + '_graphics.txt')
        if not Path(graphics_file).exists():
            raise FileNotFound(f'Graphics file {graphics_file} not found!')
        pos = open(graphics_file)
        d = json.loads(pos.read())
        conv_dict = {}
        for key, items in d.items():
            # if unique, extract the terminal modifier for later re-addition
            if self.unique:
                pattern = re.search(r'(-[0-9]+)', key)
                key = re.sub(r'(-[0-9]+)', '', key)
            try:
                conv_list = self.conversion[key]
            except KeyError:
                pass

                if self.unique:
                    # if self.unique is True, we need to add the terminal modifier back
                    conv_list = [conv + pattern.group() for conv in conv_list]
                for conv in conv_list:
                    conv_dict[conv.replace(self.prefix, '')] = items


        for key, items in d.items():
            if not key.startswith(self.species):
                conv_dict.update({key: items})
        prefix = 'up' if self.uniprot else 'ncbi-geneid'
        with open(self.wd / f'{prefix}_{Path(graphics_file).name}', 'w') as outfile:
            outfile.write(json.dumps(conv_dict))
        typer.echo(typer.style(f'Conversion of {Path(graphics_file).name} complete!', fg=typer.colors.GREEN, bold=True))



    def _process_dataframe(self, df):
        if self.unique:
            # Extract the terminal modifiers and create a new column
            # This enables the re-addition of the modifiers at the
            # end of the function.
            df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
            # Remove the terminal modifier so that the IDs map properly
            # to the KEGG API call
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)
        # Map to convert KEGG IDs to target IDs. Note lists are returned
        # for some conversions.
        df['entry1_conv'] = df['entry1'].map(self.conversion)
        df['entry2_conv'] = df['entry2'].map(self.conversion)
        # Fills nans with entries from original columns
        df['entry1'] = df['entry1_conv'].fillna(df['entry1'])
        df['entry2'] = df['entry2_conv'].fillna(df['entry2'])
        # Drop the extra column as it's all now in entry1/2 columns
        df = df.drop(['entry1_conv', 'entry2_conv'], axis=1)
        # Joins the comma seperated list to a string to avoid iterrating
        if self.uniprot and self.unique:
            df['entry1'] = [','.join(map(str, l)) for l in df['entry1']]
            df['entry2'] = [','.join(map(str, l)) for l in df['entry2']]
            # Quick seperation back into a list while avoiding cpd:, path:, and undefined
            # Also removes up: modifier to avoid breaking explode
            df['entry1'] = df['entry1'].apply(lambda x: re.findall(r'[a-zA-z0-9]+', x.replace('up:', '')) if x.startswith('up:') else [x])
            df['entry2'] = df['entry2'].apply(lambda x: re.findall(r'[a-zA-z0-9]+', x.replace('up:', '')) if x.startswith('up:') else [x])
            # Individualize each entry from a list
        df = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
        df['entry1'] = df['entry1'].str.replace(self.prefix, '')
        df['entry2'] = df['entry2'].str.replace(self.prefix, '')
        if self.unique:
            df['entry1'] = df['entry1'] + df['match1']
            df['entry2'] = df['entry2'] + df['match2']
            df = df.drop(['match1', 'match2'], axis=1)
        # Finally, remove all rows with 'hsa:' since this will create misleading files
        # Also clash with the graphics file since it won't include 'hsa:' for the for loop
        df = df[~df['entry1'].astype(str).str.startswith(self.species)]
        df = df[~df['entry2'].astype(str).str.startswith(self.species)]
        return df

    def convert_file(self):
        file = Path(self.input_data)
        df = pd.read_csv(file, delimiter='\t')
        if self.uniprot:
            typer.echo(f'Now converting {file.name} to UniProt IDs...')
            df_out = self._process_dataframe(df)
            df_out.to_csv(self.wd / 'up_{}'.format(file.name), sep='\t', index=False)
            if self.graphics != None:
                typer.echo(f'Graphics file given! Now converting {Path(self.graphics).name} to UniProt IDs...')
                self._process_graphics()
        else:
            typer.echo(f'Now converting {file.name} to NCBI IDs...')
            df_out = self._process_dataframe(df)
            df_out.to_csv(self.wd / 'ncbi_{}'.format(file.name), sep='\t', index=False)
            if self.graphics != None:
                typer.echo(f'Graphics file given! Now converting {Path(self.graphics).name} to NCBI IDs...')
                self._process_graphics()

        # print work done
        typer.echo(typer.style(f'Conversion of {file.name} complete!', fg=typer.colors.GREEN, bold=True))


def genes_convert(species, input_data, wd: Path, graphics=None,
                  uniprot: bool = False, unique: bool = False, verbose: bool = False):
    '''
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of genes that can be used in graph analysis.
    '''
    if Path(input_data).is_dir():
        for file in Path(input_data).glob('*.tsv'):
            try:
                converter = Converter(species, file, wd=wd, graphics=graphics,
                                      unique=unique, uniprot=uniprot,
                                      verbose=verbose)
                converter.convert_file()
            except FileNotFound as e:
                typer.echo(typer.style(e.message, fg=typer.colors.RED, bold=True))
                continue
    else:
        converter = Converter(species, input_data, wd, graphics=graphics,
                              unique=unique, uniprot=uniprot,
                              verbose=verbose)
        converter.convert_file()
