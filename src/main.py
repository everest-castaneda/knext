#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 20:49:11 2023

@author: everest_castaneda1
"""

import sys
import typer
import click
from pathlib import Path

from .knext.convert import convert_file
from .knext.convert import convert_folder
from .knext.call import kgml
from .knext.genes import genes_file
from .knext.genes import genes_folder
from .knext.mixed import mixed_file
from .knext.mixed import mixed_folder

@click.group()
def cli():
    pass

@cli.command()
@click.argument('species')
@click.option('-r', '--results', help = 'Directory in which to save output', required = False)
def get_kgml(species, results):
    """
    Acquires all KGML files for a given species. Use a KEGG species
    identification, usually a 3 to 4 letter organism code, as input. Handles
    directories first, hence make sure you have a proper species identifier or
    it may return empty, junk folders. The results flag is to save all results 
    to a directory. For more information about KEGG organism codes, visit: 
    https://www.genome.jp/kegg/catalog/org_list.html
    """
    if results is not None:
        results = Path(results)
        if results.exists() == False:
            typer.echo(f'Directory {results} does not exist or is invalid. Please input a valid directory...')
            sys.exit()
        else:
            kgml(species, results)
    else:
        wd = Path.cwd()
        results = wd / 'kgml_{}'.format(species)
        typer.echo(f'No output directory provided. All files will be saved to:\n{results}')
        results.mkdir(exist_ok = True)
        kgml(species, results)
        
@cli.command()
@click.argument('input_data')
@click.option('-r', '--results', required = False)
@click.option('-u', '--unique', default = False, is_flag = True)
@click.option('-g', '--graphics', default = False, is_flag = True)
@click.option('-n', '--names', default = False, is_flag = True)
def genes(input_data: str, results: str, unique: bool, graphics: bool, names: bool):
    """
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of genes that can be used in graph analysis. If -u/--unique flag 
    is used genes are returned with terminal modifiers to enhance network 
    visualization or analysis. If the -g/--graphics flag is used, x-y 
    coordinates of pathways are returned, which may be used as positions in 
    NetworkX\'s graph drawing commands.
    """
    if Path.exists(Path(input_data)) == False:
        typer.echo('Please input a directory of KGML files or an individual KGML file...')
        sys.exit()
    else:
        if results:
            if Path.exists(Path(results)) == False:
                typer.echo('Directory not found.\nPlease input a directory in which to save the results...')
            else:
                wd = Path(results)
                if graphics and unique and names:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd, unique = True, graphics = True, names = True)
                    else:
                        genes_folder(input_data, wd, unique = True, graphics = True, names = True)
                elif graphics and unique and not names:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd, unique = True, graphics = True)
                    else:
                        genes_folder(input_data, wd, unique = True, graphics = True)
                elif graphics and names and not unique:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd, graphics = True, names = True)
                    else:
                        genes_folder(input_data, wd, graphics = True, names = True)
                elif graphics and not unique and not names:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd, graphics = True)
                    else:
                        genes_folder(input_data, wd, graphics = True)
                elif unique and names and not graphics:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd, unique = True, names = True)
                    else:
                        genes_folder(input_data, wd, unique = True, names = True)
                elif unique and not names and not graphics:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd, unique = True)
                    else:
                        genes_folder(input_data, wd, unique = True)
                elif names and not unique and not graphics:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd, names = True)
                    else:
                        genes_folder(input_data, wd, names = True)
                else:
                    if Path(input_data).is_file():
                        genes_file(input_data, wd)
                    else:
                        genes_folder(input_data, wd)
                    
        else:
            wd = Path.cwd()
            typer.echo(f'\nNo output directory given. All resulting files or folders will be saved to current directory:\n{wd}\n')
            if graphics and unique and names:
                if Path(input_data).is_file():
                    genes_file(input_data, wd, unique = True, graphics = True, names = True)
                else:
                    genes_folder(input_data, wd, unique = True, graphics = True, names = True)
            elif graphics and unique and not names:
                if Path(input_data).is_file():
                    genes_file(input_data, wd, graphics = True, unique = True)
                else:
                    genes_folder(input_data, wd, graphics = True, unique = True)
            elif graphics and names and not unique:
                if Path(input_data).is_file():
                    genes_file(input_data, wd, graphics = True, names = True)
                else:
                    genes_folder(input_data, wd, graphics = True, names = True)
            elif graphics and not unique and not names:
                if Path(input_data).is_file():
                    genes_file(input_data, wd, graphics = True)
                else:
                    genes_folder(input_data, wd, graphics = True)
            elif unique and names and not graphics:
                if Path(input_data).is_file():
                    genes_file(input_data, wd, unique = True, names = True)
                else:
                    genes_folder(input_data, wd, unique = True, names = True)
            elif unique and not names and not graphics:
                if Path(input_data).is_file():
                    genes_file(input_data, wd, unique = True)
                else:
                    genes_folder(input_data, wd, unique = True)
            elif names and not unique and not graphics:
                if Path(input_data).is_file():
                    genes_file(input_data, wd, names = True)
                else:
                    genes_folder(input_data, wd, names = True)
            else:
                if Path(input_data).is_file():
                    genes_file(input_data, wd)
                else:
                    genes_folder(input_data, wd)

@cli.command()
@click.argument('input_data')
@click.option('-r', '--results', required = False)
@click.option('-u', '--unique', default = False, is_flag = True)
@click.option('-g', '--graphics', default = False, is_flag = True)
@click.option('-n', '--names', default = False, is_flag = True)
def mixed(input_data: str, results: str, unique: bool = False, graphics: bool = False, names: bool = False):
    """
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of mixed genes, compounds, and pathways that can be used in graph 
    analysis. If -u/--unique flag is used nodes are returned with terminal 
    modifiers to enhance network visualization or analysis. If the 
    -g/--graphics flag is used, x-y coordinates of pathways are returned, 
    which may be used as positions in NetworkX\'s graph drawing commands.
    """
    if Path.exists(Path(input_data)) == False:
        typer.echo('Please input a directory of KGML files or an individual KGML file...')
        sys.exit()
    else:
        if results:
            if Path.exists(Path(results)) == False:
                typer.echo('Directory not found.\nPlease input a directory in which to save the results...')
            else:
                wd = Path(results)
                if graphics and unique and names:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd, unique = True, graphics = True, names = True)
                    else:
                        mixed_folder(input_data, wd, unique = True, graphics = True, names = True)
                elif graphics and unique and not names:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd, graphics = True, unique = True)
                    else:
                        mixed_folder(input_data, wd, graphics = True, unique = True)
                elif graphics and names and not unique:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd, graphics = True, names = True)
                    else:
                        mixed_folder(input_data, wd, graphics = True, names = True)
                elif graphics and not unique and not names:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd, graphics = True)
                    else:
                        mixed_folder(input_data, wd, graphics = True)
                elif unique and names and not graphics:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd, unique = True, names = True)
                    else:
                        mixed_folder(input_data, wd, unique = True, names = True)
                elif unique and not names and not graphics:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd, unique = True)
                    else:
                        mixed_folder(input_data, wd, unique = True)
                elif names and not unique and not graphics:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd, names = True)
                    else:
                        mixed_folder(input_data, wd, names = True)
                else:
                    if Path(input_data).is_file():
                        mixed_file(input_data, wd)
                    else:
                        mixed_folder(input_data, wd)
        else:
            wd = Path.cwd()
            typer.echo(f'No output directory given. All resulting files or folders will be saved to current directory:\n{wd}\n')
            if graphics and unique and names:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd, unique = True, graphics = True, names = True)
                else:
                    mixed_folder(input_data, wd, unique = True, graphics = True, names = True)
            elif graphics and unique and not names:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd, graphics = True, unique = True)
                else:
                    mixed_folder(input_data, wd, graphics = True, unique = True)
            elif graphics and names and not unique:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd, graphics = True, names = True)
                else:
                    mixed_folder(input_data, wd, graphics = True, names = True)
            elif graphics and not unique and not names:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd, graphics = True)
                else:
                    mixed_folder(input_data, wd, graphics = True)
            elif unique and names and not graphics:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd, unique = True, names = True)
                else:
                    mixed_folder(input_data, wd, unique = True, names = True)
            elif unique and not names and not graphics:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd, unique = True)
                else:
                    mixed_folder(input_data, wd, unique = True)
            elif names and not unique and not graphics:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd, names = True)
                else:
                    mixed_folder(input_data, wd, names = True)
            else:
                if Path(input_data).is_file():
                    mixed_file(input_data, wd)
                else:
                    mixed_folder(input_data, wd)
        
@cli.command()
@click.argument('species')
@click.argument('input_data')
@click.option('-r', '--results', required = False)
@click.option('-u', '--unique', default = False, is_flag = True)
@click.option('-up', '--uniprot', default = False, is_flag = True)
@click.option('-g', '--graphics', required = False)
def convert(input_data, species, graphics, results: bool = False, uniprot: bool = False, unique: bool = False):
    """
    Converts a file or folder of parsed genes or mixed pathways from KEGG IDs to 
    NCBI gene IDs or UniProt IDs. Default is NCBI gene IDs unless the
    -up/--uniprot flag is used. The -g/--graphics flag converts the graphics file
    or folder into the accompanying IDs. Do not mix/match! Use the -u/--unique 
    flag only if the input is unique, presence of terminal modifiers. 
    Use the -g/--graphics flag with folders if the input for the conversion is 
    a folder and vice-versa.
    """
    if Path.exists(Path(input_data)) == False:
        typer.echo('Please input a directory of KGML files or an individual KGML file...')
        sys.exit()
    else:
        if results:
            wd = Path(results)
            if graphics and unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics, unique = True, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, graphics, unique = True, uniprot = True)
            elif graphics and not unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, graphics, uniprot = True)
            elif graphics and unique and not uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics, unique = True)
                else:
                    convert_folder(species, input_data, wd, graphics, unique = True)
            elif graphics and not unique and not uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics)
                else:
                    convert_folder(species, input_data, wd, graphics)
            elif not graphics and unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, unique = True, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, unique = True, uniprot = True)
            elif not graphics and not unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, uniprot = True)
            elif not graphics and unique and not uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, unique = True)
                else:
                    convert_folder(species, input_data, wd, unique = True)
            else:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd)
                else:
                    convert_folder(species, input_data, wd)       
        else:
            wd = Path.cwd()
            typer.echo(f'\nNo output directory given. All resulting files or folders will be saved to current directory:\n{wd}\n')
            if graphics and unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics, unique = True, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, graphics, unique = True, uniprot = True)
            elif graphics and not unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, graphics, uniprot = True)
            elif graphics and unique and not uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics, unique = True)
                else:
                    convert_folder(species, input_data, wd, graphics, unique = True)
            elif graphics and not unique and not uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, graphics)
                else:
                    convert_folder(species, input_data, wd, graphics)
            elif not graphics and unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, unique = True, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, unique = True, uniprot = True)
            elif not graphics and not unique and uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, uniprot = True)
                else:
                    convert_folder(species, input_data, wd, uniprot = True)
            elif not graphics and unique and not uniprot:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd, unique = True)
                else:
                    convert_folder(species, input_data, wd, unique = True)
            else:
                if Path(input_data).is_file():
                    convert_file(species, input_data, wd)
                else:
                    convert_folder(species, input_data, wd)



if __name__ == '__main__':
    cli()
