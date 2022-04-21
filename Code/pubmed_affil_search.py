from Bio import Entrez, Medline
import pandas as pd


#temp placeholders while being lazy about argparse
Entrez.email = 'adam.hockenberry@tempus.com'
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

def gather_medline_records(search_term, max_records):
    '''
    Uses esearch to retrieve medline formatted records.

    Inputs:
        search_term : str
        max_records : int
    Outputs:
        medline_records : a list
    '''
    handle = Entrez.esearch(db='pubmed', term='{}[ad]'.format(search_term), retmax=max_records)
    record = Entrez.read(handle)
    handle.close()
    id_list = record['IdList']
    assert len(id_list) == max_records 
    
    #Do a fetch to return the text data for each record
    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline",
                            retmode="text")
    medline_records = list(Medline.parse(handle))
    assert len(medline_records) == max_records
    return medline_records

def read_filter_file(filter_file):    
    '''
    Inputs:
        filter_file : str (path to existing file)
    Outputs:
        filtered_affils : list of strs
    '''
    try:
        with open(filter_file, 'r') as infile:
            filtered_affils = infile.read().splitlines()
            print('Found {} existing affiliation strings to remove from consideration:'.format(len(filtered_affils)))
    except FileNotFoundError:
        print('Did not find any previous file so starting afresh')
        filtered_affils = []
    return filtered_affils


def filter_results(search_term, medline_records, filtered_affils):
    '''

    '''
    valid_affils = []
    valid_ids = []
    for record in medline_records:
        affil_of_interest = 0
        for affil in record['AD']:
            if search_term in affil.lower():
                tempy = affil.lower()
                for error in filtered_affils:
                    tempy = tempy.replace(error, '***')
                if search_term in tempy:            
                    valid_affils.append(affil)
                    affil_of_interest += 1
        if affil_of_interest > 0:
            valid_ids.append(record['PMID'])     
    valid_affils = list(set(valid_affils))
    #Should probably write this list to a file rather than screen print
    #print('Affiliations that remain and will be considered valid:')
    #print()
    #for i in valid_affils:
    #    print(i)
    return valid_ids, valid_affils 

def fetch_xml_records(valid_ids):
    handle = Entrez.efetch(db="pubmed", id=valid_ids, rettype="medline",
                        retmode="xml")
    xml_records = Entrez.read(handle)
    return xml_records

def construct_df(xml_records, valid_affils):
    '''
    
    '''
    pmids = []
    titles = []
    dois = []
    journals = []
    years_journal = []
    volumes = []
    issues = []
    n_authors = []
    n_affil_authors = []
    affil_author_names = []
    for record in xml_records['PubmedArticle']:
        pmids.append(int(record['MedlineCitation']['PMID']))
        titles.append(record['MedlineCitation']['Article'].get('ArticleTitle'))
        journals.append(record['MedlineCitation']['Article']['Journal'].get('Title', ''))
        years_journal.append(record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''))
        volumes.append(record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Volume', ''))
        issues.append(record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Issue', ''))
        ###DOI is a bit more complicated
        doi = ''
        for item in record['MedlineCitation']['Article']['ELocationID']:
            if item.attributes['EIdType']=='doi':
                doi = str(item)
        dois.append(doi)
        ###Authorship gets really complicated
        total_authors = 0
        affil_authors = 0
        affil_author_list = []
        for person in record['MedlineCitation']['Article']['AuthorList']:
            affil_author = False
            ln = person.get('LastName', '')
            if ln == '':
                continue
            total_authors += 1
            name = (person['LastName']+' '+person['ForeName']).replace(' ', '_')
            affiliations = person['AffiliationInfo']
            for affiliation in affiliations:
                if affiliation['Affiliation'] in valid_affils:
                    affil_author=True
            if affil_author:
                affil_authors += 1
                affil_author_list.append(name)
        ###Finally adding
        n_authors.append(total_authors)
        n_affil_authors.append(affil_authors)
        affil_author_names.append('; '.join(affil_author_list))
    pubs_df = pd.DataFrame(zip(pmids, titles, dois, journals, years_journal, volumes, issues, n_authors, n_affil_authors, affil_author_names),\
            columns = ['PMID', 'Title', 'DOI', 'Journal', 'Year', 'Volume', 'Issue',\
            'Total_author_count', 'Affiliated_author_count', 'Affiliated_author_names'])
    return pubs_df

if __name__ == '__main__':
    #Obvs need to get argparser working here

    ###Read in a list of affiliation strings to remove from a file, or empty list if it doesn't exist
    filter_file = '{}/search_filters/{}_removal.txt'.format(infile_directory, '_'.join(search_term.split(' ')))
    filtered_affils = read_filter_file(filter_file)
    ###Establish how many records to grab
    max_records = get_max_pubmed_records(search_term)
    print('Getting {} records matching {}'.format(max_records, search_term))
    medline_records = gather_medline_records(search_term, max_records)
    print('Number of medline records found: {}'.format(len(medline_records)))
    valid_ids, valid_affils = filter_results(search_term, medline_records, filtered_affils)
    print('Total valid IDs to fetch XML for: {}'.format(len(valid_ids)))
    xml_records = fetch_xml_records(valid_ids)
    pubs_df = construct_df(xml_records, valid_affils)
    pubs_df = pubs_df.sort_values('Year', ascending=False)
    print('Shape of resulting dataframe: {}'.format(pubs_df.shape))
    pubs_df.to_csv('{}/{}.tsv'.format(outfile_directory, '_'.join(search_term.split(' '))), index=False, sep='\t')
