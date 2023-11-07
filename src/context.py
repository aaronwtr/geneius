from Bio import Entrez

import yaml
from tqdm import tqdm


class PubMedScraper:
    def __init__(self):
        self.config = self.config()
        self.email = self.config['user']['email']
        Entrez.email = self.email

    @staticmethod
    def config():
        with open('geneius/configs/pubmed.yml', 'r') as file:
            config = yaml.safe_load(file)
        return config

    @staticmethod
    def search_literature(query, num_records):
        try:
            search_results = Entrez.esearch(db="pubmed", term=query, retmax=num_records, sort="relevance")
            record = Entrez.read(search_results)
            pubmed_ids = record['IdList']
            return pubmed_ids
        except Exception as e:
            print(f"An error occurred: {str(e)}")
            return []

    def get_literature_context(self, query, num_records):
        pubmed_ids = self.search_literature(query, num_records)
        literature_records = []

        print("Extracting literature records...")
        for pubmed_id in tqdm(pubmed_ids):
            pubmed_record = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="text")
            subcontext = pubmed_record.read()
            abstract, title, doi = self.extract_paper_info(subcontext)
            if doi == "":
                doi = "DOI not found"
            literature_records.append([abstract, title, doi])
        context = self.combine_context(literature_records)
        return context

    @staticmethod
    def extract_paper_info(data):
        lines = data.split('\n')
        abstract = ""
        title = ""
        doi = ""
        inside_abstract = False
        inside_title = False
        inside_doi = False

        for line in lines:
            if line.startswith("AB  - "):
                inside_abstract = True
                abstract += line[6:]
            elif line.startswith("TI  - "):
                inside_title = True
                title += line[6:]
            elif line.startswith("LID - "):
                inside_doi = True
                doi += line[6:]
            elif inside_abstract:
                if line.startswith("  "):
                    abstract += " " + line.strip()
                else:
                    inside_abstract = False
            elif inside_title:
                if line.startswith("  "):  # Handle multiline titles
                    title += " " + line.strip()
                else:
                    inside_title = False
            elif inside_doi:
                if line.startswith("  "):
                    doi += " " + line.strip()
                else:
                    inside_doi = False
        return abstract, title, doi


    @staticmethod
    def combine_context(literature_records):
        claude_context = ""
        for record in literature_records:
            claude_context += f"[Title]: {record[1]} \n [DOI]: {record[2]} \n [Abstract]: {record[0]}\n\n"
        return claude_context
