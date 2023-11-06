![image](https://github.com/aaronwtr/Geneius-AnthropicAI-Hackaton/assets/54633647/42386179-410c-4711-98c0-3d18932fbd70)
# Geneius: Disease-gene hypothesis and validation powered by Anthropic AI's Claude
Geneius is a bioinformatical tool that streamlines scientific research by rapidly validating and generating a contextual understanding of gene-disease links based on simple user input prompts. Geneius can currently execute one of two specific tasks:

Task 1) Disease-Gene Link Validation: When presented with a hypothetical disease-gene link, Geneius searches through the scientific literature to find compelling evidence supporting this association. It then provides a molecular explanation, complete with a citation to the relevant research papers.

Task 2) Disease-Gene Link Hypothesis Generation: For users in need of insights about a specific disease, Geneius constructs a comprehensive disease context by extracting relevant information from the scientific literature. Claude searches and retrieves this literature, presenting users with a curated list of genes implicated in the disease based on scientific research. Furthermore, it elucidates the molecular mechanisms underpinning the involvement of these genes in the disease.

Geneius was developed at the Anthropic AI Hackathon in London, on the fourth and the fifth of November 2023. More information here: [https://devpost.com/software/geneius](https://devpost.com/software/geneius)

# Usage 
**1) Install the package and requirements**

Run `pip install geneius` and install the packages in the requirements.txt. 

<br>

**2) Generate API Key for Claude**

Go to [console.anthropic.com](console.anthropic.com) to request access to the Claude API. Once access has been granted, generate an API key and store this in a .env file in the Geneius folder. Make sure you format the Claude API secret as: 
`CLAUDE_SECRET={{your_secret}}`.  

<br>

**3) Set PubMed email in configs/pubmed.yml**

You need to provide an email when querying PubMed. You can add yours in the pubmed.yml file.

<br>

**4) Select a task**

Task 1 = Disease-gene validation, i.e. provide Geneius with a disease and a suspected linked gene, and find evidence in the literature for this link. 
Task 2 = Disease-gene hypothesis, i.e. provide Geneius with a disease and ask it to hypothesize which N genes might be underpinning this disease, with scientific substantiation.

Tasks can be selected as `--task=i` where i can be either 1 or 2.

<br>

**5) Set a disease of interest**

Store the disease you want to query the literature for in `--disease={str: your_disease}`.

<br>

**6) Set task-dependent flags**

The two different tasks require a different set of flags to allow the program to successfully execute. 

For task 1, set the following flags: 
- `--gene={str: your_gene}` (the gene you want to query)
- `--num_evidence={int: number_of_papers}` (the number of papers you want Geneius to return)

For task 2, you only need to set the flag:
- `--num_genes={int: number_of_genes}` (the number of genes you want Geneius to associate with the specified disease.

<br>

**7) Run Geneius from the command line**

Open a terminal window and run `geneius {FLAGS}`, where the flags are the task-dependent flags as outlined above.
