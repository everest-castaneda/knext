#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 21:31:31 2023

@author: everest
"""

import requests
from pathlib import Path
import typer

def kgml(species, results):
    """
    Handles the API call for acquiring the specific species input
    """
    KEGGorg = 'http://rest.kegg.jp/list/organism'
    KEGGlist = 'http://rest.kegg.jp/list/pathway/%s'
    KEGGget = 'http://rest.kegg.jp/get/%s/kgml'
    response = requests.get(KEGGorg).text.split('\n')
    org_list = []
    taxonomy = []
    # The -1 avoids the last element of the list, which is a dangling newline
    for i in range(len(response) - 1):
        org_list.append(response[i].split('\t')[1])
        taxonomy.append(response[i].split('\t')[2])
    d = {org_list[i]: taxonomy[i] for i in range(len(taxonomy))}
    if species not in org_list:
        typer.echo(f'Please input a species name in KEGG organism ID format.\nThese are usually {len(min(org_list, key = len))} to {len(max(org_list, key = len))} letter codes.\n--Example: "Homo sapiens" is "hsa"')
    else:
        typer.echo(f'Now acquiring all KGML files for {d[species]}...')
        response = requests.get(KEGGlist % species).text.split('\n')
        pathways = [r.split('\t')[0] for r in response if r]
        for path in pathways:
            typer.echo(f'Now acquiring pathway {path}...')
            config = Path(results / '{}.xml'.format(path))
            paths = requests.get(KEGGget % path).text
            if config.is_file():
                pass
            else:
                with open(results / '{}.xml'.format(path), 'w') as outfile:
                    outfile.write(paths)
            