# GBVTD

## Ground Based Velocity Track Display

This tarball contains three directories:
cappi/ : 'cappi' creator using intebilinear interpolation
findcenter/ : 'findcenter' simplex search code
vtd/ : 'VD' GBVTD wind retrieval code

To build, run 
     $ make -f Makefile.(compiler)
in each directory. 
You will need an input file for each
radar volume you which to process for VD, while findcenter supports multiple
volumes in the same input file. The input files are commented.

Notes about the Linux version:
The Linux port reads and writes files in cedric ASCII format, with the
header readable by 'Bigendian' machines (ie Solaris, HP). This is done so that
existing applications (ie cedric, grid2ps) can still read the files. It cannot
currently read or write binary cedric files,but this is planned for the future.
A reference on the cedric format can be found here:
http://www.mmm.ucar.edu/pdas/Postscript/appendix-D.ps

Good luck on the build and resulting analysis!

Report any bugs or comments to:
Michael Bell <mmbell@hawaii.edu> or
Wen-Chau Lee <wenchau@ucar.edu>