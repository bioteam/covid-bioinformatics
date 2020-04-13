#!/usr/bin/env python

import sys
import xmltodict
from Bio import Entrez

# 'Osborne BI'
search_str = sys.argv[1]
email = "briano@bioteam.net"

try:
    handle = Entrez.esearch(db="pubmed",
                            term=search_str + '[AU]',
                            retmax=1000,
                            email=email)
    result = Entrez.read(handle)
    for id in result['IdList']:
        handle = Entrez.efetch(
            db="pubmed", id=id, rettype="xml", retmode="text", email=email)
        pub = xmltodict.parse(handle.read())
        print(pub['PubmedArticleSet']['PubmedArticle']
              ['MedlineCitation']['Article']['Abstract'])

except (IOError) as exception:
    print(str(exception))

'''
pub['PubmedArticleSet']['PubmedArticle']['MedlineCitation'].keys()

['@Status', '@Owner', 'PMID', 'DateCompleted', 'DateRevised',
'Article', 'MedlineJournalInfo', 'CitationSubset', 'MeshHeadingList']


pub['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['Article'].keys()

['@PubModel', 'Journal', 'ArticleTitle', 'Pagination', 'Abstract',
            'AuthorList', 'Language', 'GrantList', 'PublicationTypeList']
'''
