# README #

### What is this repository for? ###

* KRETZ is a small library to enable reading a 3D Ultrasound Kretzfile and to convert between the toroidal geometry and cartesian space and visa versa 
* Version 1.0
* Two executabels are KretzWriter and KretzConverter to perform the conversion
* The first three tests found in the CMakeList.txt file descibe the use and command line options of the executables to convert between the different geometries. The source files for these executables show how the geometry conversion can be performed and integrated in ITK pipelines.

### How do I get set up? ###

**Dependencies**

`ITK, Boost, Doxygen and GraphViz (for documentation)`

**How to build**

Use cmake and make

**Programs**

*KretzFileWriter* either takes a KRETZ file and outputs the voxels in the geometry of the 3D ultrasound probe. One dimension corresponds to radial distance and the others correspond to the angles. The spacing is not isotropic and can be found in the KRETZ file. KretzFileWriter can take a cartesian geometry and a KRETZ file and convert the cartesian geometry back into the geometry specified in the KRETZ file.

  arguments: 
   - i - input KRETZ file
   - c - optional cartesian file to convert back into geometry specified in the KRETZ file 
   - o - path for the output file

*KretzConverter* takes a KRETZ file and converts it to cartesian coordinates.

  arguments: 
   - i - input KRETZ file
   - o - path for the output file
   - r - three floating point values corresponding to the voxel spacing in each direction, if not specified s must be
   - s - three integers to define the number of voxels in each direction, if not specifed r must be 
   - m - create a mask image 1 where the voxel is in the geometry of the ultrasound beam and 0 otherwise
   - n - normalise voxel values to have zero mean and unit variance
   - d - write out power doppler instead of grayscale voxel values to geometry of the grayscale acquisition

    
    

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
