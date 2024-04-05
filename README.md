# UE 300164-1 Advanced bioinformatic analysis of single-cell data

This repository holds the material for the course "300164-1 Advanced bioinformatic analysis of
single-cell data". It is held every summer term at the University of Vienna. Over time, I intend to
keep evolving this document and the material in the course to reflect the change in the students'
background, the needs of the department, and current developments in the analysis of single-cell
RNA-seq data.

## Course description

The main topics of the course are (accurate as of 2024 summer term):

1. Normalisation/variance stabilisation
2. Clustering and cell types
3. Integration/batch correction
4. (if time permits) Cross-species comparison

The main goal of the course is to introduce students to the theory behind scRNA-seq analysis, and
impress upon them the importance of understanding what their tools are doing. I try to put less
emphasis on reproducing results or copying what I do on the whiteboard/notebook and more on
the underlying principles.

## Changelog and notes

#### 2024

- more model species data:
    - Normalisation/variance stabilisation/intro to clustering done with [PBMC3k](https://scanpy.readthedocs.io/en/latest/generated/scanpy.datasets.pbmc3k.html#scanpy.datasets.pbmc3k) as always.
    - Real-world test with [Psilocybin orbitofrontal cortex](https://www.biorxiv.org/content/10.1101/2024.01.07.573163v1.full), where we tried to reproduce the authors' findings going off of the paper figures.
    - Not-exam: [healthy/diabetes pancreas](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061?query=E-MTAB-5061).
- including [TopOMetry](https://www.biorxiv.org/content/10.1101/2022.03.14.484134v2.full). The paper was inaccessible to students due to mathematical complexity. We tried to follow tutorials. Not sure it registered.
- including [scVI](https://scvi-tools.org); we didn't put too much emphasis on it since installation can be hard inside the conda env I ask for. Was mentioned as an option for integration.
- no cross-species comparison. We didn't have time; also the students came from a more mainstream biology background, so cross-species might have been too complicated.
- cell type discussion was difficult due to lack of common background. Maybe I should either run a short intro myself or have them all read the Opinions paper, so we can have a proper discussion.

#### 2023

(course inception)