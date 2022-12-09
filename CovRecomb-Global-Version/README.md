
# CovRecomb-Global-Version
To identify the putative inter-lineage recombinants among global sequenced SARS-CoV-2 genomes.


## What is the CovRecomb-Global-Version?
The Global-Version is designed for authors to update the global results of SARS-CoV-2 recombinants. Different from CovRcomb-Local-Version, it takes into account the epidemiology data of the analyzed genomes; thus, it could not only identify the possibility of recombination from the genomic information but also distinguish independent recombination events based on the global epidemiological background. In total, CovRecomb-Global-Version provides a semi-automatic pipeline for authors to identify recombinants and detect recombination events.


## Requirements
  - python>=3.6
  - biopython>=1.70


## Workflow
<img src="img/workflow.png"/>


### STEP 1: Transformation: transform the full genome to the mutational sites
- Date acquisition and filteration. 
- Exaction of mutations.
```
python fasta_clean_get_FM.py
python pos_norm.py
```

### STEP 2: Construction: construct the lineage-defining library
- Definition of feature mutations for each SARS-Cov-2 lineage.
- Cluster lineages with similar lineage-defining mutations.
- Parameter settings
```
python cluster_lineage.py
python least_number_of_FV.py
```

### STEP 3 & 4 & 5 
- Use the core algorithm of CovRecomb to identify inter-lineage recombinants.
- STEP 3: Predefinition: predefine a lineage-paired score matrix for each sample
- STEP 4: Mapping: map samples’ mutations against the predefined matrix
- STEP 5: Determination: determine the optimal lineage-paired combinations
```
python CovRecomb_pipeline.py
python confirm_epi_context.py
```

### STEP 6: Identification of independent recombination events
- Detect the independent recombination events from all the identified putative recombinants.
- Detect the linegae or variant paired patterns among the detected independent recombination events.
- Draw heatmap(s) for representing the lineage(and variant)-preference of recombination events.
```
python identify_patterns.py
python identify_breakpoints.py
Rscript paired_patterns.r
```

### STEP 7: Identification of the breakpoint distribution with the help of 3SEQ
- Detect the breakpoint distribution of the independent recombination events with the help of 3SEQ.

```
python identify_breakpoints.py 
```


## Acknowledgements
We sincerely thank the Global Intiative on Sharing All Influenza Data ([GISAID](https://www.gisaid.org/)) and all data contributors for making SARS-CoV-2 genomic sequence data available to the public and open science.