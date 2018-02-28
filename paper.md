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
  affiliation: "1, 3"
affiliations:
- name: Nuffield Department of Obstetrics and Gynaecology, University of Oxford,Level 3, Women’s Centre, John Radcliffe Hospital, Oxford
  index: 1
- name: School of Womens and Childrens Health, University of New South Wales, Randwick, NSW, Australia
  index: 2
- name: Fetal Medicine Unit, The Womens Centre, John Radcliffe Hospital Oxford
  index: 3
date: 14 February 2016
bibliography: paper.bib
---

# Summary

The Kretzfile format is used to store 3D ultrasound data from GE Voluson ultrasound scanners. The geometry used in these files is a toroidal coordinate system. 
We present ITK transformation and utilities to convert Kretzfiles to cartesian coordinates. Previous work [@SklicerHeart] has enabled the reading of kretz files and approximate coordinate transformations.
This work will enable medical imaging researchers to investigate clinically 3D ultrasound. 


# References
