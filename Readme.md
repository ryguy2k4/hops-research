# The Dominant Close-Companion Star Formation Pathway (Title TBD)
## Authors
Ryan Sponzilli<sup>1</sup>, Leslie W. Looney<sup>1</sup>, John J. Tobin<sup>2</sup>, Frankie J. Encalada<sup>1</sup>, Austen Fourkas<sup>1</sup>

<sup>1</sup>Department of Astronomy, University of Illinois, 1002 West Green St, Urbana, IL 61801, USA

<sup>2</sup>National Radio Astronomy Observatory, 520 Edgemont Rd., Charlottesville, VA 22903, USA

## Abstract
Understanding how close-companion protostars form is central to unraveling the processes that govern stellar multiplicity and very early star formation. We analyze an unprecedentedly large sample of Class 0/I close-companion <500 au protostellar systems using data from ALMA of 11 systems in Perseus and 40 systems in Orion. Close-companion protostars are thought to form either directly at these <500 au scales via disk fragmentation, or on >1000 au scales via turbulent fragmentation and then migrate to separations $\lesssim$500 au. The prevalence of either formation pathway is not well understood, yet it is crucial to a complete picture of the star formation process. We examine the distribution of relative position angles of the companion protostars to the position angles of the molecular outflow. The outflow is a useful proxy for the angular momentum vector of the system, orthogonal to the disk major axis. We use a simple model to account for a random sampling of inclinations and orbital phase in each system, finding that the observations are consistent with a distribution where the outflows are preferentially orthogonal to the companions. Based on this analysis, we suggest disk fragmentation is the dominant formation pathway for close-companion protostellar systems.

## Repository Structure
* `doc`: Contains project proposals, relevant references, and the research poster.
* `data`: Contains the source data used in the analysis, and also the processed data.
* `src`: Contains the source code for the reproducible analysis workflow with snakemake.
* `results`: Contains the results of the analysis, including tables and figures.
* `notebooks`: Contains all my scratch work and exploratory analysis. This is not part of the reproducible workflow, take everything in this folder with a grain of salt.

## Reproducible Workflow
The entire analysis is containerized with Docker. The `Dockerfile` in the root folder can be used to build the image. The workflow is implemented with snakemake and defined in the `Snakefile` in the root folder. The container can be run via the `runall.sh` script. Note that running the analysis requires having the relevant FITS images, which are not included in this repository. Edit `config.yaml` to point to the location of the FITS images on your system. The results of the analysis are saved in the `results` folder.

## Data
* `tobin2022_orion.txt` - sourced from Table 1: Orion Catalog in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
* `reynolds2024_perseus_1.txt` - sourced from Table 4 in https://iopscience.iop.org/article/10.3847/1538-4357/ad151d#apjad151dbib60 Table 4
* `tobin2018_perseus_2.txt` - sourced from Table 1 in https://iopscience.iop.org/article/10.3847/1538-4357/aae1f7
* `tobin2022_orion_pairings` - sourced from Table 3 in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
* `tobin2022_perseus_pairings` - sourced from Table 4 in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
* `notes.txt` - this is a google spreadsheet downloaded as a csv that contains my notes on the sources.
* FITS Images - the relevant FITS images are not included in this repository. For compatibility with the code, they should be organized such that each field is in a folder named after the field, and that folder contains a 12CO image, which contains either "12co" or "spw39" in the filename.