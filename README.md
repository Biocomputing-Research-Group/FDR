
# FineFDR

A Fine-grained Taxonomy-specific FDR assessment tool to control the FDR separately for peptides in different taxonomic ranks. The critical innovation of FineFDR is that a novel data integration model is developed to assess FDR accurately using meta-information from reference databases or sample-specific sequencing data, or from both data sources. The details are available in [paper].


## Setup

### Dependency

- python == 3.9

### Requirement

- 2 GB availbale RAM or above



## User Manul with A Toy Example

### Data Source

Download the mock microbial ”U” (UNEVEN) type community data set with the cell number U1 via the PRIDE repository PXD006118. We will call the toy example U1 in short.

http://ftp.pride.ebi.ac.uk/pride/data/archive/2017/05/PXD006118/

| File Name           | Last Modified    | Size |
|---------------------|------------------|------|
| U1_run3_1mM.raw     | 2021-11-05 19:32 | 690M |
| U1_run3_2mM.raw     | 2021-11-05 19:43 | 634M |
| U1_run3_5mM.raw     | 2021-11-05 19:24 | 606M |
| U1_run3_10mM.raw    | 2021-11-05 20:06 | 582M |
| U1_run3_20mM.raw    | 2021-11-05 19:56 | 573M |
| U1_run3_50mM.raw    | 2021-11-05 20:07 | 597M |
| U1_run3_100mM.raw   | 2021-11-05 20:08 | 572M |
| U1_run3_200mM.raw   | 2021-11-05 19:35 | 525M |
| U1_run3_500mM.raw   | 2021-11-05 19:34 | 500M |
| U1_run3_1000mM.raw  | 2021-11-05 20:04 | 474M |
| U1_run3_2000mM.raw  | 2021-11-05 19:42 | 475M |
| U1_run3_U1_11ug.raw | 2021-11-05 19:52 | 622M |

### Database Search

Search the seleted database using Comet release 2022.01 rev.0 (We will provide support for more search engines in the future releases). The sample results for U1 in .pin and .txt format are provided.

https://github.com/UWPR/Comet

### Post-search Processing (Optional)

Rerank the search result using Percolator v3.05 (We will provide support for more post-search tools in the future releases). The sample results for U1 in .pin format are provided.

https://github.com/percolator/percolator

### Generate a Specie Disctinary

Create a specie disctinary file as the following format:

```
Protein + /t + Species + /n
Protein + /t + Species + /n
Protein + /t + Species + /n
...
```
The specie disctinary for U1 in .pin and .txt format are provided.

### Run FineFDR

The list of input files include:

| File        | Description |
| ----------- | ----------- |
| Header      | Title       |
| Paragraph   | Text        |
