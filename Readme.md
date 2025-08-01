# The Dominant Close-Companion Star Formation Pathway (Title TBD)
## Authors
Ryan Sponzilli<sup>1</sup>, Leslie W. Looney<sup>1</sup>, John J. Tobin<sup>2</sup>, Frankie J. Encalada<sup>1</sup>

<sup>1</sup>Department of Astronomy, University of Illinois, 1002 West Green St, Urbana, IL 61801, USA

<sup>2</sup>National Radio Astronomy Observatory, 520 Edgemont Rd., Charlottesville, VA 22903, USA

## Abstract
We analyze an unprecedentedly large sample of class 0/I close-companion <500 AU systems using data from ALMA of 10 fields in Perseus and 40 fields in Orion. Close-companion protostars are thought to form either directly at these <500 AU scales via disk fragmentation, or they form on >1000 AU scales via turbulent fragmentation and then migrate to separations <500 AU. It is of interest to know which mechanism dominates. We identify and measure outflows by eye, and use outflow alignments to constrain the formation mechanism. We compare cumulative frequency distributions of ΔPA (difference between the outflow PA and binary separation PA), and the observations most closely match what we would expect to see if ΔPA is always 90°. We find that disk fragmentation is the dominant formation pathway for close-companion protostellar systems, and turbulent fragmentation followed by inward migration is a less common formation pathway.

## Repository Structure
* `doc`: Contains project proposals, relevant references, and the research poster.
* `data`: Contains the source data used in the analysis, and also the processed data.
* `src`: Contains the source code for the reproducible analysis workflow with snakemake.
* `results`: Contains the results of the analysis, including tables and figures.
* `notebooks`: Contains all my scratch work and exploratory analysis. This is not part of the reproducible workflow, take everything in this folder with a grain of salt.

## Snakemake Workflow
The worflow is implemented with snakemake and defined in the `Snakefile` in the root folder. With all the requirements installed (see `requirements.txt`), the worflow can be executed in the terminal with `snakemake runall`.

## Data
* `tobin2022_orion.txt` - sourced from Table 1: Orion Catalog in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
* `reynolds2024_perseus_1.txt` - sourced from Table 4 in https://iopscience.iop.org/article/10.3847/1538-4357/ad151d#apjad151dbib60 Table 4
* `tobin2018_perseus_2.txt` - sourced from Table 1 in https://iopscience.iop.org/article/10.3847/1538-4357/aae1f7
* `tobin2022_orion_pairings` - sourced from Table 3 in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
* `tobin2022_perseus_pairings` - sourced from Table 4 in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
* `notes.txt` - this is a google spreadsheet downloaded as a csv that contains my notes on the sources.
* FITS Images - the relevant FITS images are not included in this repository. For compatibility with the code, they should be organized such that each field is in a folder named after the field, and that folder contains a 12CO image, which contains either "12co" or "spw39" in the filename.