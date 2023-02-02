# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining all KGML files for a single species
"""

import typer
import requests
from pathlib import Path

app = typer.Typer()

@app.command()
def kgml(species_name: str = typer.Argument(None, file_okay = False, dir_okay = False, help = 'Species name in KEGG organism code format. Visit https://www.genome.jp/kegg/catalog/org_list.html for more information.')):
    KEGGorg = 'http://rest.kegg.jp/list/organism'
    KEGGlist = 'http://rest.kegg.jp/list/pathway/%s'
    KEGGget = 'http://rest.kegg.jp/get/%s/kgml'
    response = requests.get(KEGGorg).text.split('\n')
    org_list = []
    taxonomy = []
    for i in range(len(response) - 1):
        org_list.append(response[i].split('\t')[1])
        taxonomy.append(response[i].split('\t')[2])
    d = {org_list[i]: taxonomy[i] for i in range(len(taxonomy))}
    if species_name not in org_list:
        typer.echo(f'Please input a species name in KEGG organism ID format.\nThese are usually {len(min(org_list, key = len))} to {len(max(org_list, key = len))} letter codes.\n--Example: "Homo sapiens" is "hsa"')
    else:
        typer.echo(f'Now acquiring all KGML files for {d[species_name]}...')
        wd = Path.cwd()
        output_path = wd / 'kgml_{}'.format(species_name)
        output_path.mkdir(exist_ok = True)
        response = requests.get(KEGGlist % species_name).text
        pathways = [r.replace('path:', '') for r in response.split() if r.startswith('path:')]
        for path in pathways:
            typer.echo(f'Now acquiring pathway {path}...')
            config = Path(output_path / '{}.xml'.format(path))
            paths = requests.get(KEGGget % path).text
            if config.is_file():
                pass
            else:
                with open(output_path / '{}.xml'.format(path), 'w') as outfile:
                    outfile.write(paths)
            
if __name__ == '__main__':
    app()
