# README #

### What is this repository for? ###

* KRETZ is a small library to enable reading a 3D Ultrasound Kretzfile and to convert between the toroidal geometry and cartesian space and visa versa 
* Version 0.01
* Two executabels are KretzWriter and KretzConverter to perform the conversion
* The first three tests found in the CMakeList.txt file descibe the use and command line options of the executables to convert between the different geometries. The source files for these executables show how the geometry conversion can be performed and integrated in ITK pipelines.

### How do I get set up? ###

**Dependencies**

`ITK, Boost, Doxygen and GraphViz (for documentation)`

**How to build**

Use cmake and make

**How to run tests**

make test

**How to create docs**

make doc

## Community Guidelines

Please use the [issue tracker](https://github.com/plooney/kretz/issues) to report any problems with the software. If you want to contribute to the development of KRETZ, please send a pull request.

## Contributors

* Padraig Looney, University of Oxford

## License

KRETZ is provided under the terms of the [GNU Lesser General Public License, version 3](https://www.gnu.org/licenses/lgpl-3.0.en.html). 

## Acknowledgements

* This project has received funding from was supported by the Eunice Kennedy Shriver National Institute of Child Health and Human Development (NICHD) Human Placenta Project of the National Institutes of Health under award number UO1-HD087209.
