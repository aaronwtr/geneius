import requests
from bs4 import BeautifulSoup
import os

class PubMedAbstractScraper:
    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(self, gene_list):
        self.gene_list = gene_list

    def search_pubmed_for_gene(self, gene_name):
        params = {
            'db': 'pubmed',
            'term': gene_name,
            'retmode': 'json',
            'retmax': 250  
        }
        response = requests.get(self.ESEARCH_URL, params=params)
        
        if response.status_code == 200:
            pmids = response.json()['esearchresult']['idlist']
            return pmids
        else:
            print(f"Error searching PubMed: {response.status_code}")
            return []

    def fetch_details_of_pmids(self, pmids):
        ids = ','.join(pmids)
        params = {
            'db': 'pubmed',
            'id': ids,
            'retmode': 'xml',
        }
        response = requests.get(self.EFETCH_URL, params=params)
        
        if response.status_code == 200:
            return BeautifulSoup(response.content, 'xml')
        else:
            print(f"Error fetching details: {response.status_code}")
            return None

    def extract_abstracts(self, xml_data):
        abstracts = {}
        articles = xml_data.find_all('PubmedArticle')
        for article in articles:
            pmid = article.find('PMID').text
            abstract_text = article.find('AbstractText')
            if abstract_text:
                abstracts[pmid] = abstract_text.get_text()
        return abstracts

    def score_abstracts(self, abstracts, gene_name):
        scored_abstracts = []
        for pmid, abstract in abstracts.items():
            gene_count = abstract.upper().count(gene_name.upper())
            scored_abstracts.append((gene_count, pmid, abstract))
        scored_abstracts.sort(reverse=True, key=lambda x: x[0])
        return scored_abstracts[:10]

    def save_abstracts_to_file(self, scored_abstracts, gene_name, filename='top_abstracts.txt'):
        directory = 'data'
        if not os.path.exists(directory):
            os.makedirs(directory)
        filepath = f"{directory}/{gene_name}_{filename}"
        with open(filepath, 'w') as file:
            for score, pmid, abstract in scored_abstracts:
                file.write(f"PMID: {pmid}\nScore: {score}\nAbstract: {abstract}\n\n")

    def run(self):
        for gene_name in self.gene_list:
            print(f"Processing gene: {gene_name}")
            pmids = self.search_pubmed_for_gene(gene_name)
            if pmids:
                xml_data = self.fetch_details_of_pmids(pmids)
                if xml_data:
                    abstracts = self.extract_abstracts(xml_data)
                    scored_abstracts = self.score_abstracts(abstracts, gene_name)
                    self.save_abstracts_to_file(scored_abstracts, gene_name)
                    for score, pmid, _ in scored_abstracts:
                        print(f"Gene: {gene_name}, PMID: {pmid}, Score: {score}")

# Usage
gene_list = ['BRCA1', 'TP53', 'EGFR']  # Replace with the list of genes of interest
scraper = PubMedAbstractScraper(gene_list)
scraper.run()

# import requests
# from bs4 import BeautifulSoup

# class PubMedAbstractScraper:
#     ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
#     EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

#     def __init__(self, gene_name):
#         self.gene_name = gene_name

#     def search_pubmed_for_gene(self):
#         params = {
#             'db': 'pubmed',
#             'term': self.gene_name,
#             'retmode': 'json',
#             'retmax': 250  # Increased to have a larger pool for scoring
#         }
#         response = requests.get(self.ESEARCH_URL, params=params)
        
#         if response.status_code == 200:
#             pmids = response.json()['esearchresult']['idlist']
#             return pmids
#         else:
#             print(f"Error searching PubMed: {response.status_code}")
#             return []

#     def fetch_details_of_pmids(self, pmids):
#         ids = ','.join(pmids)
#         params = {
#             'db': 'pubmed',
#             'id': ids,
#             'retmode': 'xml',
#         }
#         response = requests.get(self.EFETCH_URL, params=params)
        
#         if response.status_code == 200:
#             return BeautifulSoup(response.content, 'xml')
#         else:
#             print(f"Error fetching details: {response.status_code}")
#             return None

#     def extract_abstracts(self, xml_data):
#         abstracts = {}
#         articles = xml_data.find_all('PubmedArticle')
#         for article in articles:
#             pmid = article.find('PMID').text
#             abstract_text = article.find('AbstractText')
#             if abstract_text:
#                 abstracts[pmid] = abstract_text.get_text()
#         return abstracts

#     def score_abstracts(self, abstracts):
#         scored_abstracts = []
#         for pmid, abstract in abstracts.items():
#             gene_count = abstract.upper().count(self.gene_name.upper())
#             scored_abstracts.append((gene_count, pmid, abstract))
#         scored_abstracts.sort(reverse=True)
#         return scored_abstracts[:10]

#     def save_abstracts_to_file(self, scored_abstracts, filename='top_abstracts.txt'):
#         with open(filename, 'w') as file:
#             for score, pmid, abstract in scored_abstracts:
#                 file.write(f"PMID: {pmid}\nScore: {score}\nAbstract: {abstract}\n\n")

#     def run(self):
#         pmids = self.search_pubmed_for_gene()
#         if pmids:
#             xml_data = self.fetch_details_of_pmids(pmids)
#             if xml_data:
#                 abstracts = self.extract_abstracts(xml_data)
#                 scored_abstracts = self.score_abstracts(abstracts)
#                 self.save_abstracts_to_file(scored_abstracts)
#                 for score, pmid, _ in scored_abstracts:
#                     print(f"PMID: {pmid}, Score: {score}")

# # Usage
# gene_name = 'BRCA1'  # Replace with the gene of interest
# scraper = PubMedAbstractScraper(gene_name)
# scraper.run()
