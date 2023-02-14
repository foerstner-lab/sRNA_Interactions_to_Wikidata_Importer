# sRNA interactions to Wikidata Importer
[![License](https://img.shields.io/pypi/l/reademption.svg)](https://github.com/foerstner-lab/sRNA_Interactions_to_Wikidata_Importer)
[![DOI](https://zenodo.org/badge/18210971.svg)](https://zenodo.org/record/7638542)
## Description
This tool can be used to import interactions between small RNAs and genes to Wikidata. This tool was developed as a part of [InteractOA](https://interactoa.toolforge.org/) project. Currently, it only support importing data from [RegulonDB](https://regulondb.ccg.unam.mx/).
It automates the importing process
## Future development
Initially, this tool was developed to show the application of small RNA data modelling in Wikidata, and was developed for parsing RegulonDB data only. Now, we are intending to further extend it and make generic.

## Prerequisites
1. Install dependencies:
   - ```pip install wikidataintegrator pywikibot bcbio-gff pandas```
2. Modify user-config.py
    - Go to line 43
    - Add your Wikidata's username
## Usage
1. Run ```python run.py```
2. The running script will ask for further inputs
