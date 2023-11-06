from context import PubMedScraper
from models import claude
from anthropic import HUMAN_PROMPT, AI_PROMPT

import time
import math
import argparse


def main():
    parser = argparse.ArgumentParser(description="Geneius: A tool for biomedical literature search and extraction.")
    parser.add_argument("--task", type=int, choices=[1, 2], required=True, help="Task number (1 or 2)")
    parser.add_argument("--disease", type=str, required=True, help="Disease name")
    parser.add_argument("--gene", type=str, help="Gene name (only for Task 1)")
    parser.add_argument("--num_evidence", type=int, help="Number of evidence (only for Task 1)")
    parser.add_argument("--num_genes", type=int, help="Number of genes (only for Task 2)")
    args = parser.parse_args()

    start_time = time.time()

    pms = PubMedScraper()
    disease = args.disease

    _claude = None
    prompt = None
    response = None
    num_records = None

    if args.task == 1:
        if args.gene is None or args.num_evidence is None:
            print("For Task 1, both --gene and --num_evidence are required.")
            return
        gene = args.gene
        num_evidence = args.num_evidence

        context, num_records = pms.get_literature_context(disease)
        _claude = claude.Claude()
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

    elif args.task == 2:
        if args.num_genes is None:
            print("For Task 2, --num_genes is required.")
            return
        num_genes = args.num_genes

        context, num_records = pms.get_literature_context(disease)
        _claude = claude.Claude()
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
    response = _claude.create_completion(prompt)

    print(response)

    print(f"Collected and parsed through {num_records} scientific papers in: "
          f"{(math.floor((time.time() - start_time) / 60))} minutes and {math.floor((time.time() - start_time) % 60)} "
          f"seconds.")
