# Young Binary Formation Mechanisms: Outflow Analysis
## Authors
Ryan Sponzilli<sup>1</sup>, Leslie W. Looney<sup>1</sup>, John J. Tobin<sup>2</sup>, Frankie J. Encalada<sup>1</sup>

<sup>1</sup>Department of Astronomy, University of Illinois, 1002 West Green St, Urbana, IL 61801, USA

<sup>2</sup>National Radio Astronomy Observatory, 520 Edgemont Rd., Charlottesville, VA 22903, USA

## Overview

### Introduction
We investigate the origins of stellar multiplicity by examining very young binary and multiple-star systems. Two formation mechanisms are proposed to account for multiplicity: disk fragmentation (e.g., Kratter et al. 2010) and turbulent fragmentation (e.g., Oﬀner et al. 2010). As a young star forms, it accretes material from its collapsing envelope, and some is channeled along the poles, driving an outflow of material. Disk fragmentation predicts that outflows will be orthogonal to the disk plane and turbulent fragmentation predicts random orientations.

### Observations
This project uses data from the Atacama Large Millimeter Array (ALMA), a radio interferometer in Chile. We have 40 sources from the Orion molecular cloud and 11 from the Perseus molecular cloud. Each source recorded data for multiple molecular emission lines; we analyze the 12CO molecular line data since 12CO is an excellent tracer of protostellar outflows, and after H2, 12CO is the most abundant molecule beyond Earth. Of the 51 close-multiple systems, 36 had visible outflow structures. 

### Methods
The outflow position angles were measured manually. Outflows are typically v-shaped structures and may be visible on one or both sides of the protostar. We measured the position angle of the center of the outflow and averaged them if there were two sides. We calculate ΔPA, the difference between the measured outflow position angle and the observed separation position angle. The separation between the binaries lies parallel to their orbital plane, so ΔPA is a proxy for outflow alignment. 

### Results
We construct a cumulative frequency distribution to compare multiple possible models. Each model is constructed with a simple simulation. Our observed distribution of ΔPA most closely matches a model that predicts all outflows to be orthogonal. The distribution includes angles less than 90° which is an effect of the inclination of the system and its projection on the sky. This result indicates that disk fragmentation is the dominant formation mechanism within our sample.

## Repository Structure
* `doc`: Contains project proposal, relevant references, and the research poster.
* `data`: Contains the source data used in the analysis, and also the processed data.
* `src`: Contains the source code for the reproducible analysis workflow with snakemake.
* `results`: Contains the results of the analysis, including figures of the outflow measurements, and plots illustrating the measured data.
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