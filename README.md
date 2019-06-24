# Intramolecular cross-linking

Intramolecular cross-links are created within a single macromolecular chain MD model on the basis of distance between cross-linking sites. Cross-linking sites are randomly chosen in advance and a pair of cross-linking sites falling within the capture radius is cross-linked. Intramolecular cross-links are allowed to form at a distance larger than bond equilibrium length and are relaxed to the bond equilibrium length through multistep equilibration.

## Getting Started

### crosslink.py
Generates and runs LAMMPs input scripts.

### identify.py
Chooses cross-linking sites.

## Running
```
python crosslink.py 1000 10
```
crosslink.py requires two variables. The first one is the degree of polymerization of the given linear model and the second one is the target degree of cross-linking.

# Reference
@article{10.1039/C7SM00360A,	
author = {Bae, Suwon and Galant, Or and Diesendruck, Charles E. and Silberstein, Meredith N.},	
title = {Tailoring single chain polymer nanoparticle thermo-mechanical behavior by cross-link density},	
journal = {Soft Matter},	
volume = {13},	
issue = {15},	
pages = {2808-2816},	
year = {2017},
doi = {10.1039/C7SM00360A},
URL = {
http://dx.doi.org/10.1039/C7SM00360A
}

}
