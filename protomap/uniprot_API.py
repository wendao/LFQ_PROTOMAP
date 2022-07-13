#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file: uniprot_API.py
@comment: 
@created: 2022/04/14 20:35:05
@auther: Zhu, Yinyue
@version: 1.0
'''
import os
import pandas as pd
from urllib import parse,request
from io import StringIO

# for more parameters, refer to https://www.uniprot.org/help/api_queries
params = {
'from': 'ID', # identifier, refer to https://www.uniprot.org/help/api_idmapping
'to': 'SWISSPROT',
'format': 'tab', # output format
#'query': '', # identifiers split by white characters('\w')
'columns':'id,entry_name,protein_names,genes', # refer to https://www.uniprot.org/help/uniprotkb_column_names
#'taxon':'9606', # refer to https://www.uniprot.org/taxonomy/
}


def id_mapping(query,params=params):
    '''
    Mapping identifier list or file to another identifier format.
    '''
    url = 'https://www.uniprot.org/uploadlists/'
    try:
        if os.path.isfile(query):
            file=open(query)
            params['query']=file.read()
    except:
        params['query']=' '.join(query)
    data = parse.urlencode(params)
    data = data.encode('utf-8')
    req = request.Request(url, data)
    with request.urlopen(req) as f:
       response = f.read()
    return pd.read_csv(StringIO(response.decode('utf-8')),sep='\t')
    #print(response.decode('utf-8'))