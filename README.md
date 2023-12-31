![image](https://github.com/aaronwtr/Geneius-AnthropicAI-Hackaton/assets/54633647/42386179-410c-4711-98c0-3d18932fbd70)
# Geneius: Disease-gene hypothesis and validation powered by Anthropic AI's Claude
Geneius is a bioinformatics command-line tool that streamlines scientific research by rapidly validating and generating a contextual understanding of gene-disease links based on simple user input prompts. Geneius can currently execute one of two specific tasks:

Task 1) Disease-Gene Link Validation: When presented with a hypothetical disease-gene link, Geneius searches through the scientific literature to find compelling evidence supporting this association. It then provides a molecular explanation, complete with a citation to the relevant research papers.

Task 2) Disease-Gene Link Hypothesis Generation: For users in need of insights about a specific disease, Geneius constructs a comprehensive disease context by extracting relevant information from the scientific literature. Claude searches and retrieves this literature, presenting users with a curated list of genes implicated in the disease based on scientific research. Furthermore, it elucidates the molecular mechanisms underpinning the involvement of these genes in the disease.

Tasks can be selected by specifying args. For task 1, you will need to specify a `--gene`, for task 2 you will need to specify `--num_genes`.

Geneius was developed at the Anthropic AI Hackathon in London, on the fourth and the fifth of November 2023. More information here: [https://devpost.com/software/geneius](https://devpost.com/software/geneius)

# Usage 
**1) Generate API Key for Claude**

Go to [console.anthropic.com](console.anthropic.com) to request access to the Claude API. Once access has been granted, generate an API key and store this in a .env file in the Geneius folder. Make sure you format the Claude API secret as: 
`CLAUDE_SECRET={{your_secret}}`.  

<br>

**2) Install the package and requirements**

Run `pip install geneius` and install the packages in the requirements.txt. 

<br>

**3) Set PubMed email in configs/pubmed.yml**

You need to provide an email when querying PubMed. You can add yours in the pubmed.yml file.

<br>

**4) Put in your Anthropic API key (REQUIRED)**

Add the Anthropic API key in the `--api_key` flag. 

<br>

**5) Set a disease flag (REQUIRED)**

Store the disease you want to query the literature for in `--disease={str: your_disease}`.

<br>

**6) Set a num_records flag (REQUIRED)**

Decide how many scientific papers Geneius will look through. Claude's 100k token context window permits a maximum of about 600 papers. `--num_records={int: num_records}`.

<br>

**7) Set task-dependent flags**

The two different tasks require a different set of flags to allow the program to successfully execute. 

For task 1, set the following flags: 
- `--gene={str: your_gene}` (the gene you want to query)

For task 2, you only need to set the flag:
- `--num_genes={int: number_of_genes}` (the number of genes you want Geneius to associate with the specified disease.

<br>

**8) Run Geneius from the command line**

Open a terminal window and run `geneius {FLAGS}`, where the flags are the task-dependent flags as outlined above.
