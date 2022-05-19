from Bio import Entrez, Medline
import pandas as pd
import os

from dotenv import load_dotenv
load_dotenv()

Entrez.email = os.environ.get('PUBMED_EMAIL') 
search_term = 'tempus'
#search_term = 'guardant'
infile_directory = '../Data'
outfile_directory = '../Results'



def get_max_pubmed_records(search_termm, db='pubmed'):
    '''
    Owing to some intricacies of pubmed, subsequent sesarches have a retmax argument so we need to specify
    the expected number of hits first to ensure we retrieve all records.
    
    Inputs:
        search_term : str

    Outputs:
        max_records : int

    '''
    max_records = 0
    handle = Entrez.egquery(term='{}[ad]'.format(search_term))
    record = Entrez.read(handle)
    handle.close()
    for row in record['eGQueryResult']:
        if row['DbName'] == db:
            max_records = int(row['Count'])
    return max_records

def gather_ids_all_hits(search_term, max_records):
    '''
    Uses esearch to retrieve medline formatted records.

    Inputs:
        search_term : str
        max_records : int
    Outputs:
        id_list : a list
    '''
    handle = Entrez.esearch(db='pubmed', term='{}[ad]'.format(search_term), retmax=max_records)
    record = Entrez.read(handle)
    handle.close()
    id_list = record['IdList']
    return id_list 

def fetch_xml_records(valid_ids):
    '''
    Bulk fetching of xml records from id list using pubmed

    Inputs:
        valid_ids : list of ids (ints)

    Outputs:
        xml_records : list of xml records

    '''
    handle = Entrez.efetch(db="pubmed", id=valid_ids, rettype="medline",
                        retmode="xml")
    xml_records = Entrez.read(handle)
    handle.close()
    return xml_records

def add_individual_record(record, pubs_df, next_index):
    '''

    '''
    ###Easy stuff
    pubs_df.at[next_index, 'PMID'] = record['MedlineCitation']['PMID']
    pubs_df.at[next_index, 'Title'] = record['MedlineCitation']['Article'].get('ArticleTitle')
    pubs_df.at[next_index, 'Journal'] = record['MedlineCitation']['Article']['Journal'].get('Title', '')
    pubs_df.at[next_index, 'Year'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', '')
    pubs_df.at[next_index, 'Volume'] = record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Volume', '')
    pubs_df.at[next_index, 'Issue'] = record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Issue', '')
    
    ###DOI is a bit more complicated
    doi = ''
    for item in record['MedlineCitation']['Article']['ELocationID']:
        if item.attributes['EIdType']=='doi':
            doi = str(item)
    pubs_df.at[next_index, 'DOI'] = doi

    ###Authorship information gets unfortunately complicated
    total_authors = 0
    affil_authors = 0
    affil_author_list = []
    affil_list = []
    for person in record['MedlineCitation']['Article']['AuthorList']:
        affil_author = False
        ln = person.get('LastName', '')
        if ln == '':
            continue
        total_authors += 1
        name = (person['LastName']+' '+person['ForeName']).replace(' ', '_')
        affiliations = person['AffiliationInfo']
        for affiliation in affiliations:
            if search_term in affiliation.lower():
                affil_author=True
                affil_hit = affiliation
        if affil_author:
            affil_authors += 1
            affil_author_list.append(name)
            affil_list.append(affil_hit)
    pubs_df.at[next_index, 'Total_author_count'] = n_authors
    pubs_df.at[next_index, 'Affiliated_author_count'] = n_affil_authors
    pubs_df.at[next_index, 'Affiliated_author_names'] = '; '.join(affil_author_list)
    pubs_df.at[next_index, 'Affiliated_author_affiliations'] = '; '.join(affil_list)

    return pubs_df

def extend_df(xml_records, pubs_df):
    '''
    
    Inputs:

    Outputs:


    '''
    new_pmids = []
    for record in xml_records['PubmedArticle']:
        if int(record['MedlineCitation']['PMID']) in pubs_df['PMID']:
            continue
        next_index = pubs_df.shape[0]
        new_pmids.append(record['MedlineCitation']['PMID'])
        pubs_df = add_individual_record(record, pubs_df, next_index)
    return pubs_df, new_pmids

def create_blank_tsv(outfile_path):
    pubs_df = pd.DataFrame(columns = ['PMID', 'Title', 'DOI', 'Journal', 'Year', 'Volume', 'Issue',\
            'Total_author_count', 'Affiliated_author_count', 'Affiliated_author_names', 'Affilated_author_affiliations', 'Is_valid_hit'])
    pubs_df.to_csv(outfile_path, sep='\t')
    return


if __name__ == '__main__':
    #Obvs need to get argparser working here
    

    

    ################
    search_term_clean = '_'.join(search_term.split(' ')) 

    
    ###Grab the available dataframe or start a blank one
    try:
        pubs_df = pd.read_csv('{}/{}.tsv'.format(outfile_directory, search_term_clean), sep='\t')
    except FileNotFoundError:
        create_blank_tsv('{}/{}.tsv'.format(outfile_directory, search_term_clean))
        pubs_df = pd.read_csv('{}/{}.tsv'.format(outfile_directory, search_term_clean), sep='\t')

    ###Establish how many records to fetch
    max_records = get_max_pubmed_records(search_term)
    ###Get the relevant IDs
    id_list = gather_ids_all_hits(search_term, max_records)
    ###Fetch the XMLs
    xml_records = fetch_xml_records(id_list)
    
    ###Extend the existing dataframe
    pubs_df, new_pmids = extend_df(xml_records, pubs_df)
    ###Re-sort
    pubs_df = pubs_df.sort_values('Year', ascending=False)
    ###And save
    pubs_df.to_csv('{}/{}.tsv'.format(outfile_directory, search_term_clean, index=True, sep='\t')
