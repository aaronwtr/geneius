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
        with open('../configs/pubmed.yml', 'r') as file:
            config = yaml.safe_load(file)
        return config

    @staticmethod
    def search_literature(query, retmax=50):
        try:
            search_results = Entrez.esearch(db="pubmed", term=query, retmax=retmax, sort="relevance")
            record = Entrez.read(search_results)
            pubmed_ids = record['IdList']
            return pubmed_ids
        except Exception as e:
            print(f"An error occurred: {str(e)}")
            return []

    def get_literature_context(self, query):
        pubmed_ids = self.search_literature(query)
        literature_records = []

        for pubmed_id in tqdm(pubmed_ids):
            pubmed_record = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="text")
            subcontext = pubmed_record.read()
            abstract, title = self.extract_abstract_and_title(subcontext)
            literature_records.append([abstract, title])
        context = self.combine_context(literature_records)
        return context

    @staticmethod
    def extract_abstract_and_title(data):
        lines = data.split('\n')
        abstract = ""
        title = ""
        inside_abstract = False
        inside_title = False

        for line in lines:
            if line.startswith("AB  - "):
                inside_abstract = True
                abstract += line[6:]
            elif line.startswith("TI  - "):
                inside_title = True
                title += line[6:]
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
        return abstract, title


    @staticmethod
    def combine_context(literature_records):
        claude_context = ""
        for i, record in enumerate(literature_records):
            claude_context += f"{i + 1} [Title]: {record[1]} \n [Abstract]: {record[0]}\n\n"
        return claude_context


if __name__ == '__main__':
    pms = PubMedScraper()

    claude_context = pms.get_literature_context("TP53")
    print(len(claude_context))