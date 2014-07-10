synthetic_datasets_nspheres
===========================

An n-spheres based synthetic data generator for supervised classification for Matlab.

# Introduction

This code implements the synthetic data generator described in:

J. Sánchez-Monedero, P.A. Gutiérrez, M. Pérez-Ortiz, and C. Hervás-Martínez. *An n-spheres based synthetic data generator for supervised classification*. In Ignacio Rojas, Gonzalo Joya, and Joan
Gabestany, editors, Advances in Computational Intelligence. 12th International Work-Conference on Artificial Neural Networks, IWANN 2013, volume 7902 of Lecture Notes in Computer Science,
pages 613–621. Springer, 2013. URL: (http://dx.doi.org/10.1007/978-3-642-38679-4_62)

## Abstract of the paper

Synthetic datasets can be useful in a variety of situations, specifically when new machine learning models and training algorithms are developed or when trying to seek the weaknesses of an specific method. In contrast to real-world data, synthetic datasets provide a controlled environment for analysing concrete critic points such as outlier tolerance, data dimensionality influence and class imbalance, among others. In this paper, a framework for synthetic data generation is developed with special attention to pattern order in the space, data dimensionality, class overlapping and data multimodality. Variables such as position, width and overlapping of data distributions in the n-dimensional space are controlled by considering them as n-spheres. The method is tested in the context of ordinal regression, a paradigm of classification where there is an order arrangement between categories. The contribution of the paper is the full control over data topology and over a set of relevant statistical properties of the data. 

# Example datasets for 1-3 dimensions:

## K = 1, sigma = 0.125
<img src="http://www.uco.es/ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-1-sig-0.125000-modes5.png" mce_src="/grupos/ayrna/../../ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-1-sig-0.125000-modes5.png" alt="Example image representing syntetic datasets" width="500" />

## K = 1, sigma = 0.125

<img src="http://www.uco.es/ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-1-sig-0.125000-modes7" mce_src="/grupos/ayrna/../../ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-1-sig-0.125000-modes7" alt="Example image representing syntetic datasets" width="500" />

## K = 2, sigma = 0.250

<img src="http://www.uco.es/ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-2-sig-0.250000-modes3.png" mce_src="/grupos/ayrna/../../ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-2-sig-0.250000-modes3.png" alt="Example image representing syntetic datasets" width="500" />

## K = 3, sigma = 0.125 (without displaying the n-spheres)

<img src="http://www.uco.es/ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.125000-modes1.png" mce_src="/grupos/ayrna/../../ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.125000-modes1.png" alt="Example image representing syntetic datasets" width="500" />

## K = 3, sigma = 0.125 (including n-spheres)

<img src="http://www.uco.es/ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.125000-modes1-sd.png" mce_src="/grupos/ayrna/../../ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.125000-modes1-sd.png" alt="Example image representing syntetic datasets" width="500" />

## K = 3, sigma = 0.250 </b><b>(without displaying the n-spheres)

<img src="http://www.uco.es/ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.250000-modes7.png" mce_src="/grupos/ayrna/../../ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.250000-modes7.png" alt="Example image representing syntetic datasets" width="500" />

## K = 3, sigma = 0.250 </b><br /><b>(including n-spheres)

<img src="http://www.uco.es/ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.250000-modes7-sd.png" mce_src="/grupos/ayrna/../../ayrna/sourcecode/iwann2013-syntheticdatagenerator/Synthetic-Gaussian-K-3-sig-0.250000-modes7-sd.png" alt="Example image representing syntetic datasets" width="500" />

# Citation and license 	

If you use this software please cite it propertly with the following bitex entry:

> @INPROCEEDINGS{SanchezMonedero2013iwann,
>  author = {J. S\'anchez-Monedero and P.A. Guti\'errez and M. P\'erez-Ortiz and C. Herv\'as-Mart\'inez},
>  title = {An {\itshape n}-Spheres Based Synthetic Data Generator for Supervised Classification},
>  booktitle = {Advances in Computational Intelligence. 12th International Work-Conference on Artificial Neural Networks, IWANN 2013},
> year = {2013},
> editor = {Ignacio Rojas and Gonzalo Joya and Joan Gabestany},
> volume = {7902},
> series = {Lecture Notes in Computer Science},
> pages = {613--621},
> publisher = {Springer},
> isbn = {978-3-642-38678-7},
> location = {Heidelberg}
> }

Please, send bugs and feedback to jsanchezm at uco dot es. All the code is GPLv3 licenced with the exception of the external tools included (plot2svg and export_fig are included).


