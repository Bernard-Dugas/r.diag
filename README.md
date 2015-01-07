
 R.DIAG - A GEM model and analysis Diagnostics Toolkit
 ======

 An extensive toolkit available to users wanting to manipulate 2D
 and 3D meteorological data, as produced either by the Environment
 Canada (EC) GEM forecast model or the EC CCCma General Circulation
 Model (GCM). The CMC/RPN GEM file formats and two flavors of the
 GCM's binary sequential format are supported.

 The toolkit is built in the R.DIAG executable binary. The toolkit
 commands can be grouped under several themes or section.

 The available sections are :

 a) Display operations
 b) File/label/record manipulations 
 c) General computations
 d) Manipulations on gridded data
 e) Manipulations on Fourrier or Spherical harmonics data
 f) Time-series manipulations

 In addition to the toolkit, a conversion tool for either the CCCma
 or CMC/RPN file formats to/from the NetCDF v3 file format is included.
 Input NetCDF files should closely conform to the CF v1.4 Metadata
 convention. The converter will otherwise (at best) choke on them.
 Depending on this executable's name, which should be either cdf2ccc
 or cdf2rpn, the default EC file format read/written by the executable
 will be either the CCCma of CMC/RPN formats. These two executable
 are in general automatically generated as hard-links.

 To generate the toolkit executable, the RPN/CMC development environment
 has to be installed and active (see mfvalin/rmnlib on github.com). As
 well, the Vgrid Descriptors package used in the GEM v4+ model also has
 to be available. Furthermore, to generate the two cdf2xxx executables,
 the old NetCDF v3.6 and UdUnits v1.2 library packages also have to
 available.

 The Linux version of the toolkit also uses the DDFUN90 package as
 produced by  David H. Bailey of the NERSC, Lawrence Berkeley Lab.
 The 2005-03-11 version of this package is included here.
 
 R.DIAG is copyrighted (C) 1990-2010 by the "Division de Recherche
 en Prevision Numerique" of Environment Canada. This code is free
 software; you can redistribute it and/or modify it under the terms
 of the GNU Lesser General Public License as published by the Free
 Software Foundation, version 2.1 of the License.

 Contact : Dugas.Bernard@uqam.ca
 Last revision : January 2015

