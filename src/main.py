from context import PubMedScraper
from models import claude
from anthropic import HUMAN_PROMPT, AI_PROMPT


if __name__ == '__main__':
    pms = PubMedScraper()
    disease = "colorectal cancer"
    gene = "APC"

    context = pms.get_literature_context(disease)

    claude = claude.Claude()
    prompt = f"{HUMAN_PROMPT}Imagine you are a researcher going through the literature to extract valuable information " \
             f"for implicating involvement of gene {gene} in disease {disease}. In order to effectively retrieve " \
             f"information, I will provide with with context from scientific literature. You can use your internal data" \
             f"and this context to formulate a response. If you are uncertain, do not speculate. ONLY return a list of " \
             f"DOI links confirming the connection of the disease and the gene, if there are any. If there is no DOI found," \
             f"return the paper title. The output should look like: 'DOI link 1' \n 'DOI link 2' ... Do not provide any" \
             f" additional explanatory sentences, just the paper title's and/or DOI's if available. \n\n" \
             f"[Context]: {context} {AI_PROMPT}"
    response = claude.create_completion(prompt)
    print(response)
