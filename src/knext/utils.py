import re
from collections import defaultdict
import urllib.request as request

import numpy as np


def conv_dict(root):
    '''
    Parse "entry" elements in the KEGG API XML file for the given root and returns a dictionary
    with the entry id as the key and the entry name as the value.
    '''
    # This dictionary is the default option
    # Only compounds are unique
    entry_id, entry_name, entry_type = _parse_entries(root)

    unique_compound = []
    for i in range(0, len(entry_id)):
        c = [name + '-' + entry_id[i] if name.startswith('cpd:') or name == 'undefined' else name for name in entry_name[i].split()]
        unique_compound.append(' '.join(c))

    conversion_dictionary = dict(zip(entry_id, unique_compound))
    return conversion_dictionary

def conv_dict_unique(root):
    # This dictionary is the unique version
    # Every item is unique to reveal subgraphs
    entry_id, entry_name, entry_type = _parse_entries(root)

    unique_name = []
    for i in range(0, len(entry_id)):
        e = [name + '-' + entry_id[i] for name in entry_name[i].split()]
        unique_name.append(' '.join(e))

    conversion_dictionary = dict(zip(entry_id, unique_name))
    return conversion_dictionary

def graphics_dict(root):
    '''
    Parses the graphics in the KEGG API XML file for the given root.
    Returns list of tuple of the entry id and the x and y coordinates.
    '''
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
            try:
                # Adds to dictionary the end entry, which is the written out name
                dd[n] = re.sub('^ ', '', list(s)[0].split(';')[1])
            except IndexError:
                # Some genes only have a name and no description
                dd[n] = split_response[0].split('\t')[1]
        # Only obtains compounds
        elif n.startswith('cpd:'):
            # Remove terminal modifiers, which are always added to compounds
            # unless a non-unique mixed pathway is chosen
            n4url = re.sub(r'-[0-9]+', '', n)
            # Uses find to get gene info since other api tools give error
            url = 'https://rest.kegg.jp/find/compound/%s'
            response = request.urlopen(url % n4url).read().decode('utf-8')
            subbed_response = re.sub(r'%s\t' % n4url, '', response)
            try:
                # Find only the query compound if given back several accessions that are similar
                split_response = re.sub('^ ', '', subbed_response.strip('\n').split(';')[1])
            except IndexError:
                # Some compounds only have one name
                split_response = subbed_response.strip('\n')
            # Adds to dictionary the end entry, which is the written out name
            dd[n] = split_response
        elif n.startswith('path:'):
            n4url1 = re.sub(r'-[0-9]+', '', n)
            n4url2 = re.sub(r'path:{}'.format(organism), '', n4url1)
            url = 'https://rest.kegg.jp/find/pathway/%s'
            response = request.urlopen(url % n4url2).read().decode('utf-8').strip('\n').split('\t')
            try:
                dd[n] = response[1]
            except IndexError:
                # One pathway has no metadata and gives an error if line not included
                dd[n] = np.nan
        else:
            dd[n] = np.nan
    return dd

def _parse_entries(root):
    '''
    Parses the entries in the KEGG API XML file for the given root.
    Returns the entry id, name, and type.
    '''
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

    return entry_id, entry_name, entry_type