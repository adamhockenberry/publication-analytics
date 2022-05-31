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

def gather_pmids(search_term, max_records):
    '''
    Uses esearch to retrieve pubmed IDs (PMID) 

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

def add_individual_record(record, pubs_df):
    '''
    Add rows to existing dataframe. This is not the ideal way to grow a dataframe but is computationally tractable
    enough to be inefficient and slow. Dataframe construction is far faster using growing lists and converting
    them into a df but it's clunky.

    Inputs:
        record : xml record (returned from pubmed search result)
        pubs_df : existing dataframe object

    Outputs:
        pubs_df : updated dataframe

    '''
    ###Easy stuff
    pmid = record['MedlineCitation']['PMID']
    pubs_df.at[pmid, 'Title'] = record['MedlineCitation']['Article'].get('ArticleTitle')
    pubs_df.at[pmid, 'Journal'] = record['MedlineCitation']['Article']['Journal'].get('Title', '')
    pubs_df.at[pmid, 'Year'] = int(record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''))
    pubs_df.at[pmid, 'Volume'] = record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Volume', '')
    pubs_df.at[pmid, 'Issue'] = record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Issue', '')
    
    ###DOI is a bit more complicated
    doi = ''
    for item in record['MedlineCitation']['Article']['ELocationID']:
        if item.attributes['EIdType']=='doi':
            doi = str(item)
    pubs_df.at[pmid, 'DOI'] = doi

    ###Authorship information gets unfortunately complicated
    total_authors = 0
    affil_authors = 0
    affil_author_list = []
    affil_list = []
    for person in record['MedlineCitation']['Article']['AuthorList']:
        affil_author = False
        ln = person.get('LastName', '')
        if ln == '': ###skip and ignore individual if I can't find a last name
            continue
        total_authors += 1
        name = (person['LastName']+' '+person['ForeName']).replace(' ', '_')
        affiliations = person['AffiliationInfo']
        for affiliation in affiliations:
            if search_term in affiliation['Affiliation'].lower():
                affil_author=True
                affil_hit = affiliation['Affiliation'].lower()
        if affil_author:
            affil_authors += 1
            affil_author_list.append(name)
            affil_list.append(affil_hit)
    pubs_df.at[pmid, 'Total_author_count'] = total_authors
    pubs_df.at[pmid, 'Affiliated_author_count'] = affil_authors
    pubs_df.at[pmid, 'Affiliated_author_names'] = '; '.join(affil_author_list)
    pubs_df.at[pmid, 'Affiliated_author_affiliations'] = '; '.join(affil_list)

    return pubs_df

def extend_df(xml_records, pubs_df):
    '''
    Iterate to add list of records to existing df. Will not update existing records at all

    Inputs:
        xml_records : list of xml
        pubs_df : dataframe
    
    Outputs:
        pubs_df : extended dataframe
        new_pmids : list of ints to track additions

    '''
    new_pmids = []
    for record in xml_records['PubmedArticle']:
        if int(record['MedlineCitation']['PMID']) in pubs_df.index:
            continue ###Skip any entries that already exist
        next_index = pubs_df.shape[0]
        new_pmids.append(record['MedlineCitation']['PMID'])
        pubs_df = add_individual_record(record, pubs_df)
    return pubs_df, new_pmids

def create_blank_tsv(outfile_path):
    pubs_df = pd.DataFrame(columns = ['PMID', 'Title', 'DOI', 'Journal', 'Year', 'Volume', 'Issue',\
            'Total_author_count', 'Affiliated_author_count', 'Affiliated_author_names', 'Affilated_author_affiliations', 'Is_valid_hit'])
    pubs_df.set_index('PMID', inplace=True)
    pubs_df.to_csv(outfile_path, index=True, sep='\t')
    return


if __name__ == '__main__':
    #Obvs need to get argparser working here
    

    

    ################
    search_term_clean = '_'.join(search_term.split(' ')) 

    
    ###Grab the available dataframe OR start a blank one
    try:
        pubs_df = pd.read_csv('{}/{}.tsv'.format(outfile_directory, search_term_clean), index_col='PMID', sep='\t')
    except FileNotFoundError:
        create_blank_tsv('{}/{}.tsv'.format(outfile_directory, search_term_clean))
        pubs_df = pd.read_csv('{}/{}.tsv'.format(outfile_directory, search_term_clean), index_col='PMID', sep='\t')

    ###Establish how many records to fetch
    max_records = get_max_pubmed_records(search_term)
    ###Get the relevant IDs
    id_list = gather_pmids(search_term, max_records)
    ###Fetch the XMLs
    xml_records = fetch_xml_records(id_list)
    
    ###Extend the existing dataframe
    pubs_df, new_pmids = extend_df(xml_records, pubs_df)
    ###Re-sort
    pubs_df['Year'] = pubs_df['Year'].astype(pd.Int64Dtype())
    pubs_df = pubs_df.sort_values('Year', ascending=False, na_position='last')
    ###And save
    pubs_df.to_csv('{}/{}.tsv'.format(outfile_directory, search_term_clean), index=True, sep='\t')
