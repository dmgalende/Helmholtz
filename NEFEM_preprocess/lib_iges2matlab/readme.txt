


			COLLECTION OF IGES FUNCTIONS FOR MATLAB
			=======================================
			    	      version 1.3
			    	      ===========



This is a version of the free collection of IGES functions for MATLAB. It can be downloaded
for free at Mathworks file exchange community site. This version can not handle all IGES entity types
but I hope this collections of functions will be of great help to do calculations on entities for
you in your IGES file. To add more entities in iges2matlab and the other functions read and
understand the IGES Specification version5x.pdf found at

http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=12087&objectType=file

in the ZIP-file. If you would like to share your upgraded version with other, please send it to
me, per.bergstrom@ltu.se . 

Important to know is that the free NURBS toolbox must be installed or saved in the same directory
as this collection of files. The NURBS toolbox can be downloaded at

http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=312&objectType=file

nrbeval from NURBS toolbox is used in this collection.
Compiled versions for Windows, Linux and Solaris is submited in zip-file. Other
users must download the correct version first.



				MAIN FUNCTIONS
				==============


IGES2MATLAB
-----------

Function for converting IGES file format to Matlab data format. Can not handle all IGES entity types.



PLOTIGES
--------

Plots lines, curves, points and surfaces from IGES file.



TRANSFORMIGES
-------------

Transform the parameter data in IGES file with a rotation/reflection matrix and a translation vector.



PROJIGES
--------

Returns points of projections on surfaces from an IGES-file.



PROJPARTIGES
------------

Returns points of projections on part of surfaces from an IGES-file.



				SUBFUNCTIONS
				============


RETSRFCRVPNT
------------

A function which returns values from surfaces, curves and points. No complete documentation is
given.



BSPEVAL.DLL/MEX
---------------

Function from the NURBS tolbox written by Mark Spink and can be downloaded for free at

http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=312&objectType=file




For more documentation about the functions above, see the help for each function in Matlab.

If you like this file collection and have use of it, please let me know that.




/ Per Bergström

per.bergstrom@ltu.se




 