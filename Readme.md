# Protostellar Outflows Shed Light on the Dominant Close Companion Star Formation Pathways
## Authors
Ryan Sponzilli<sup>1</sup>, Leslie W. Looney<sup>1</sup>, John J. Tobin<sup>2</sup>, Frankie J. Encalada<sup>1</sup>, Austen Fourkas<sup>1</sup>, Hector Arce<sup>3</sup>, Erin Cox<sup>4,5</sup>, James Di Francesco<sup>6,7</sup>, Nicole Karnath<sup>8</sup>, Zhi-Yun Li<sup>9</sup>, Nadia Murillo<sup>10</sup>, Stella Offner<sup>11</sup>, Sarah Sadavoy<sup>12</sup>, Rajeeb Sharma<sup>13</sup>

<sup>1</sup> Department of Astronomy, University of Illinois, 1002 West Green St, Urbana, IL 61801, USA
<sup>2</sup> National Radio Astronomy Observatory, 520 Edgemont Rd., Charlottesville, VA 22903, USA
<sup>3</sup> Department of Astronomy, Yale University, P.O. Box 208101, New Haven, CT 06520, USA
<sup>4</sup> NSF-Simons AI Institute for the Sky (SkAI), 172 E. Chestnut St., Chicago, IL 60611, USA
<sup>5</sup> Center for Interdisciplinary Exploration and Research in Astrophysics (CIERA), 1800 Sherman Avenue, Evanston, IL 60201, USA
<sup>6</sup> NRC Herzberg Astronomy and Astrophysics, 5071 West Saanich Road, Victoria, BC V9E 2E7, Canada
<sup>7</sup> Department of Physics and Astronomy, University of Victoria, Victoria, BC V8P 5C2, Canada
<sup>8</sup> Space Science Institute, 4765 Walnut Street, Suite B, Boulder, CO 80301, USA
<sup>9</sup> Department of Astronomy, University of Virginia, P.O. Box 400325, 530 McCormick Road, Charlottesville, VA 22904-4325, USA
<sup>10</sup> Instituto de Astronomía, Universidad Nacional Autónoma de México, AP106, Ensenada, CP 22830, B.C., Mexico
<sup>11</sup> Department of Astronomy, University of Texas at Austin, TX 78712, USA
<sup>12</sup> Department of Physics and Astronomy, York University, Toronto, Ontario M3J 1P3, Canada
<sup>13</sup> Niels Bohr Institute, University of Copenhagen, Jagtvej 155A, 2200 Copenhagen N., Denmark

## Abstract
Understanding the formation pathway for close-companion protostars is central to unraveling the processes that govern stellar multiplicity and very early star formation. We analyze a large sample of 51 Class 0/I close-companion protostellar systems where 38 of these systems show detectable outflows, yielding 42 measured outflows used in our analysis. We use ALMA observations of 11 systems in Perseus and 40 systems in Orion. These companions formed either directly at these small scales (<~ 500 au separations) via disk fragmentation or at larger scales (>1000 au separations) via turbulent fragmentation followed by inward migration. Because of differences in formation mechanism, the former is expected to have preferentially aligned disks and outflows and the latter is expected to have non-preferentially aligned disks and outflows. The relative prevalence of these formation pathways remains uncertain, yet it is critical to forming a comprehensive picture of star formation.  We examine the distribution of position angles of companion protostars relative to the position angles of their molecular outflows. The outflow, as traced by 12CO (J=2-->1), is a useful proxy for the angular momentum of the system, expected to be orthogonal to the binary orbital plane. We use a simple model to account for a random sampling of inclination and orbital phase in each system, finding that the observations are consistent with a distribution where the outflows are preferentially orthogonal to the companions. Based on this analysis, we suggest disk fragmentation is the dominant formation pathway for close-companion protostellar systems.

## Repository Structure
* `runall.sh` : Runs the workflow in a docker container if all input data is present.
* `Dockerfile` : Specifies the Docker image.
* `Snakefile` : Specifies the workflow tasks.
* `requirements.txt` : Specifies the python libraries needed to run the workflow.
* `config.yaml` : Specifies the location of the FITS files.
* `notebooks`: Contains scratch work and exploratory analysis.
* `doc`: Contains project proposals.
* `data`
    * `input` : Folder contains tables sourced from references with source positions, separations, distances, and protostellar classes for the systems in this analysis.
        * `tobin2022_orion.txt` : Table 1 in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
        * `reynolds2024_perseus_1.txt` : Table 4 in https://iopscience.iop.org/article/10.3847/1538-4357/ad151d#apjad151dbib60
        * `tobin2018_perseus_2.txt` : Table 1 in https://iopscience.iop.org/article/10.3847/1538-4357/aae1f7
        * `tobin2022_orion_pairings` : Table 3 in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
        * `tobin2022_perseus_pairings` : Table 4 in https://iopscience.iop.org/article/10.3847/1538-4357/ac36d2
        * `notes.txt` :  Contains my notes and measurements of each system
        * `offsets.csv` : Contains manually determined RA and DEC offsets for each source.
    * `output` : Folder contains organized csv tables generated in `prepare_tables.py` which are used in the analysis and presented in the paper.
* `scripts`
    * `_create_figs.py` : Defines various functions that assist with creating the figures produced in the analysis.
    * `make_all_m0_outflow_maps.py` : Generates `results/m0_outflow_maps.pdf`. 
    * `make_all_m8_maps.py` : Generates `results/all_m8_maps.pdf`.
    * `make_figs_for_paper.py` : Generates the files in `results/figs_for_paper`
    * `master_reference_doc.py` : Generates `results/master_reference.pdf`.
    * `prepare_tables.py` : Reorganizes the input data into the output data and generates latex tables presented in the paper.
    * `stat_test-combined.py` : Performs the statistical analysis described in Section 3.3 of the paper and generates the files in `results/stat_test`.
* `results`
    * `all_m8_maps.pdf` : Contains all the M8 maps, one per page, for convenience.
    * `m0_outflow_maps.pdf` : Contains all the M0 outflow maps, one per page, for convenience.
    * `master_reference.pdf` : Contains all the continuum maps, M8 maps, and M0 maps for each system, one system per page, for convenience.
    * `figs_for_paper` : Contains 15 pdf files of figures which contain M0 or M8 maps which are presented in the paper.
    * `stat_test` : Contains the DeltaPA histogram figure presented in the paper, the DeltaPA cumulative frequency distribution presented in the paper, and a figure plotting p-value vs % orthogonal distribution which is not presented in the paper.
    * `tables` : Contains the 4 latex formatted tables presented in the paper.

## Reproducible Workflow
The entire analysis is containerized with Docker. The `Dockerfile` in the root folder can be used to build the image. The workflow is implemented with snakemake and defined in the `Snakefile` in the root folder. The container can be run via the `runall.sh` script. Note that running the analysis requires having the relevant FITS images. Edit `config.yaml` to point to the location of the FITS images on your system. The results of the analysis are saved in the `results` folder.

## FITS Images
The relevant FITS images are not included in this repository. For compatibility with the code, they should be organized such that each field is in a folder named after the field, and that folder contains a 12CO image, which contains either "12co" or "spw39" in the filename.
