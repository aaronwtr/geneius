from context import PubMedScraper
from models import claude
from anthropic import HUMAN_PROMPT, AI_PROMPT

import time
import math


if __name__ == '__main__':
    task = 1

    # TASK 1: Given a disease and a gene, generate a disease context from the literature and use it to generate a
    # response to elucidate the molecular mechanism of the gene's involvement in the disease.

    if task == 1:
        start_time = time.time()
        pms = PubMedScraper()
        disease = "colorectal cancer"
        gene = "APC"
        num_evidence = 5

        context, num_records = pms.get_literature_context(disease)

        claude = claude.Claude()
        prompt = f"{HUMAN_PROMPT} Imagine you are an expert researcher going through the literature to extract " \
                 f"{num_evidence} pieces of evidence implicating molecular involvement of gene {gene} in disease " \
                 f" {disease}. I want you to explain the molecular mechanism of the gene's involvement in " \
                 f"the disease based on the scientific context I am providing you. In order to " \
                 f"effectively retrieve information, I will provide you with context from scientific literature. You " \
                 f"can use your internal data and this context to formulate a response. If you are uncertain, do not " \
                 f"speculate. Restrict yourself to returning information confirming the connection of the disease  " \
                 f"and the gene, if there are any. Strictly return only papers that have a DOI available. Your  " \
                 f"response should look like <response>[Title]: 'paper title'\n [DOI]:'doi'\n [Explanation]: This " \
                 f"paper suggests [gene] is linked to [disease] [reason]</response> Take care to complete all " \
                 f"fields of your response entirely. \n\n  <context>{context}</context> {AI_PROMPT}"
        response = claude.create_completion(prompt)

        print(response)

        print(f"Collected and parsed through {num_records} scientific papers in: "
              f"{(math.floor((time.time()-start_time)/60))} minutes and {math.floor((time.time()-start_time)% 60)} "
              f"seconds.")

    # TASK 2: Given a disease, find N genes that are implicated in the disease and generate a response to elucidate
    # the molecular mechanism of the genes' involvement in the disease.

    if task == 2:
        start_time = time.time()
        pms = PubMedScraper()
        disease = "colorectal cancer"
        num_genes = 5

        context, num_records = pms.get_literature_context(disease)

        claude = claude.Claude()
        prompt = f"{HUMAN_PROMPT} Imagine you are an expert researcher going through the literature to find " \
                 f"{num_genes} genes that are involved in {disease}, and corresponding evidence implicating  " \
                 f"molecular involvement of the genes in disease {disease}. I want you to explain " \
                 f"the molecular mechanism of the gene's involvement in " \
                 f"the disease based on the scientific context I am providing you. In order to " \
                 f"effectively retrieve information, I will provide you with context from scientific literature. You " \
                 f"can use your internal data and this context to formulate a response. If you are uncertain, do not " \
                 f"speculate. Restrict yourself to returning information confirming the connection of the " \
                 f"disease and the gene, if there are any. Strictly only restrict to papers that have a DOI available" \
                 f". Your response  should look like <response>[Genes]: [Gene 1, Gene 2, ... Gene N] \n [Title]: " \
                 f"'paper title'\n [DOI]:'doi'\n [Explanation]: This paper suggests [gene] is linked to [disease] " \
                 f"[reason]</response> Take care to complete all fields of your response entirely. \n\n" \
                 f"<context>{context}</context> {AI_PROMPT}"
        response = claude.create_completion(prompt)

        print(response)

        print(f"Collected and parsed through {num_records} scientific papers in: "
              f"{(math.floor((time.time()-start_time)/60))} minutes and {math.floor((time.time()-start_time)% 60)} "
              f"seconds.")
