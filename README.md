
# FineFDR

A fast and accurate Fine-grained Taxonomy-specific FDR assessment tool to control the FDR separately for peptides in different taxonomic ranks. 
The details are available in [paper].


## Setup

### Dependency

- python == 3.9

### Requirement

- 2 GB availbale RAM or above

## User Manul with A Toy Example

### Run FineFDR in 2 steps

#### Step 1 - Customize the config file
1. FineFDR config file for FDR control: An example config file is provided as fconfig.cfg.

```
# Decoy Pre-fix
decoy_prefix=Rev_
# FDR at PSM level
psm_fdr=0.01
# FDR at Peptide level
pep_fdr=0.01
# FDR at Protein level
protein_fdr=0.01
# The path to comet search result folder
comet_txt_dir=./toy_example/hgut_txt/
# Mehod: "C" for Comet (e-value), "CP" for Comet + Percolator (Percolator score)
method=C
# Path to the taxonomy database file
sp_dic=./toy_example/humanGutMarine_ProToOTU.tsv
# The path to the percolator search result target pin file
percolator_target=./toy_example/target_hgut.tsv
# The path to the percolator search result decoy pin file
percolator_decoy=./toy_example/decoy_hgut.tsv
# The path for the output result
output_dir=./output_result
```
2. Sipros-Ensemble config file for protein assembly: An example config file is provided as SiprosConfig.cfg.

https://github.com/guo-xuan/Sipros-Ensemble/tree/master/configs

#### Step 2 - Command line

```
python fdr_main.py [your config file]
```
For example,

```
python fdr_main.py fconfig.cfg
```


### Data Source of the Toy Example

Download the human gut microbial complex via the PRIDE repository PXD006118. We will call the toy example HG in short.

Or, download its database result directly (https://github.com/Biocomputing-Research-Group/FDR/blob/main/toy_example.7z). And, extract files.

### Database Search

Search the seleted database using Comet release 2022.01 rev.0 (We will provide support for more search engines in the future releases). 

The sample results for HG in .pin and .txt format are provided.

https://github.com/UWPR/Comet

### Post-search Processing (Optional)

Rerank the search result using Percolator v3.05 (We will provide support for more post-search tools in the future releases). 

The sample results for HG in .pin format are provided.

https://github.com/percolator/percolator

### Generate a taxonomy database file

Create a taxonomy database file file in the following format:

```
Protein + /t + Species + /n
Protein + /t + Species + /n
Protein + /t + Species + /n
...
```
The taxonomy database for HG are provided.

### IO

The list of input files:

| Files        | Description |
| ----------- | ----------- |
| *.txt     | Comet search results in txt format      |
| *.pin   | Percolator re-ranked result in pin format        |
| *.tsv     |Taxonomy database file      |
| *.fasta   | Database        |

The details of output files:

https://github.com/guo-xuan/Sipros-Ensemble/blob/master/OUTPUT.md

## Feedback

If you have any feedback, please reach out to us at Xuan.Guo@unt.edu.


## License

[GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html)

