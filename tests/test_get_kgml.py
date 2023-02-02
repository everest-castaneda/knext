#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 20:31:40 2023

@author: everest_castaneda1
"""

import pytest
import requests
from pathlib import Path
import xml.etree.ElementTree as ET

def test_get_kgml():
    KEGGget = 'http://rest.kegg.jp/get/hsa00232/kgml'
    response = requests.get(KEGGget).text
    root = ET.fromstring(response)
    
    tree_data = ET.parse('data/hsa00232.xml')
    root_data = tree_data.getroot()
    
    assert root_data.get('name') == root.get('name')
    