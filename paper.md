---
title: '3D ultrasound file reading and coordinate transformations'
tags:
- ITK
- 3D Ultrasound
authors:
- name: Pádraig T. Looney
  orcid: 0000-0002-0764-5413
  affiliation: 1 
- name: Gordon N. Stevenson
  orcid: 0000-0002-2809-8084
  affiliation: 2
- name: Sally L. Collins
  orcid: 0000-0002-0648-7433
  affiliation: 1
affiliations:
- name: University of Oxford
  index: 1
- name: University of New South Wales
  index: 2
date: 14 February 2016
bibliography: paper.bib
---

# Summary

The Kretzfile format is used to store 3D ultrasound data from GE Voluson ultrasound scanners. The geometry used in these files is a toroidal coordinate system. 
We present ITK transformation and utilities to convert Kretzfiles to cartesian coordinates. Previous work [@SklicerHeart] has enabled the reading of kretz files and approximate coordinate transformations.
This work will enable medical imaging researchers to investigate clinically 3D ultraound. 


# References
