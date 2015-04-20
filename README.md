
 R.DIAG - A GEM model and analysis Diagnostics Toolkit
 ======

 An extensive toolkit available to users wanting to manipulate 2D
 and 3D meteorological data, as produced either by the Environment
 Canada (EC) GEM forecast model or the EC CCCma General Circulation
 Model (GCM). The CMC/RPN GEM file formats and two flavors of the
 GCM's binary sequential format are supported. This toolkit was
 originally derived from an early 1990's unix port of the then
 CRAY version of the CCCma diagnostic toolkit. Additional code
 produced at UQAM for diagnostics of regional climate data in
 the late 1990's is also included in the toolbox.

 The GEM toolkit is built in the R.DIAG executable binary. The
 toolkit commands can be grouped under several themes or section.

 The available sections are :

 a) Display operations
 b) File/label/record manipulations 
 c) General computations
 d) Manipulations on gridded data
 e) Manipulations on Fourrier or Spherical harmonics data
 f) Time-series manipulations

 Once installed, documentation .txt files and a few .html files can
 be found in the /man/pdoc sub-directory. Each R.DIAG commands is thus
 documented. In particular, the /man/pdoc/index.html file holds a
 desciption of the packages and a few of the basic arguments and
 environment variables that it responds to. The .txt files are
 only available in English while the .html file are in French
 at this time.

 In addition to the toolkit, a conversion tool for either the CCCma
 or CMC/RPN file formats to/from the NetCDF v3 file format is included.
 Input NetCDF files should closely conform to the CF v1.4 Metadata
 convention. The converter will otherwise (at best) choke on them.
 Depending on this executable's name, which should be either cdf2ccc
 or cdf2rpn, the default EC file format read/written by the executable
 will be either the CCCma or CMC/RPN formats. These two executable
 are in general automatically generated as hard-links. The first
 version of the converter was created by the Ouranos Consortium
 from 2003 to 2006.

 To generate the toolkit executable, the RPN/CMC development environment
 has to be installed and active (see mfvalin/rmnlib on github.com). As
 well, the Vgrid Descriptors package used in the GEM v4+ model also has
 to be available. Furthermore, to generate the two cdf2xxx executables,
 the old NetCDF v3.6 and UdUnits v1.2 library packages also have to
 available. Most of the code available here is written in FORTRAN.

 Some versions of the toolkit may also use the DDFUN90 package as
 produced by  David H. Bailey of the NERSC, Lawrence Berkeley Lab.
 The 2005-03-11 version of this package is included here. This is
 the case for the PGI versions as these compilers (at least as of
 their version 14xx) do not provide for quad-precision real arithmetic.
 On the other hand, since the Intel, AIX/xlf and GFORTRAN compilers do
 support this, version generated with them will not require the
 DDFUN90 package. The file src/lssub/gaussg.F90 may need to be
 modified according to your environment's specifications.
 
 R.DIAG is copyrighted (C) 1990-2010 by the "Division de Recherche
 en Prevision Numerique" of Environment Canada. This code is free
 software; you can redistribute it and/or modify it under the terms
 of the GNU Lesser General Public License as published by the Free
 Software Foundation, version 2.1 of the License.

 Contact : Dugas.Bernard@uqam.ca
 Last revision : April 2015

