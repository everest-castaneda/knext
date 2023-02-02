#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 09:55:10 2023

@author: everest_castaneda1
"""

import xml.etree.ElementTree as ET
import json
import pandas as pd
import networkx as nx

def test_parser():
    tree = ET.parse('data/hsa00232.xml')
    root = tree.getroot()
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

    edgelist = []
    for line in d:
        x = line[0].replace("{","").replace("}","").replace('"entry1":',"")\
            .replace('"entry2":',"").replace('"type":',"").replace('"name":',"")\
            .replace('"value":',"").replace('"','').replace(',',"\t").replace(' ', '')
        y = line[1].replace("{","").replace("}","").replace('"name":',"")\
            .replace('"',"").replace('value:',"").replace(",",'\t').replace(' ', '')
        edgelist.append(x+"\t"+y)

    df = pd.DataFrame(edgelist)
    df = df[0].str.split("\t", expand=True).rename({0: 'entry1',1: 'entry2',
                                                  2: 'types', 3:'name',
                                                  4: 'value'}, axis='columns')
    G = nx.from_pandas_edgelist(df, source = 'entry1', target = 'entry2')
    assert G.edges