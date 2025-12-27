# Maximal Empty Spheres Contouring

Basic code for SDF contouring, as described in our paper [A Polyhedral Construction of Empty Spheres in Discrete Distance Fields](https://dl.acm.org/doi/10.1145/3721238.3730748).

This repository contains two versions of the reconstruction:
- (recommended) a cpp implementation based on CGAL,
- (deprecated) the (updated) python implementation used in our paper.

For instructions on how to build the cpp version see the Readme in the `cgal/` folder.

## Example Notebooks
We provide two demo notebooks, make sure to run `jupyter notebook` from the root directory (here).

Both notebooks require to build the cpp version first.

- `contouring_example.ipynb`: a minimal example of the application
- `contouring_comparison.ipynb`: a comparison with other contouring approaches
.

# Citation
If you use this code in your research, please cite 

```
@inproceedings{10.1145/3721238.3730748,
author = {Kohlbrenner, Maximilian and Alexa, Marc},
title = {A Polyhedral Construction of Empty Spheres in Discrete Distance Fields},
year = {2025},
isbn = {9798400715402},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3721238.3730748},
doi = {10.1145/3721238.3730748},
abstract = {Lie sphere geometry provides a unified representation of points, oriented spheres and hyperplanes in Euclidean d-space as the subset of lines in ( mathbb {R}^{d+3} ) that are contained in a certain quadric. The natural scalar product in this construction is zero if two elements are in oriented contact. We show how the sign of this product can be used to decide if spheres are disjoint. This allows us to model the space of spheres that are not intersecting a given union of spheres as the intersection of half-spaces (and the quadric). The maximal spheres are on the boundary of this set and can be computed by first constructing the intersection of half-spaces, which is a convex hull problem, and then intersecting edges of the hull against the quadric, which are the roots of a univariate quadratic. We demonstrate the method at the example of contouring a discrete signed distance field: every sample of the signed distance field represents an empty spheres and the zero-level contour has to be disjoint from the union of these spheres. Maximal spheres outside the empty spheres provide samples on the zero-level contour. The quality of this sample set is comparable to existing methods relying on optimization, while being deterministic and faster in practice.},
booktitle = {Proceedings of the Special Interest Group on Computer Graphics and Interactive Techniques Conference Conference Papers},
articleno = {22},
numpages = {10},
keywords = {Lie sphere geometry, polyhedra, discrete distance fields, maximal empty spheres},
location = {
},
series = {SIGGRAPH Conference Papers '25}
}
```
