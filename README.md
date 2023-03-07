MAKE_SURVEY
===========

The `make_survey` code is meant to take a periodic box of objects 
(galaxy mocks, halos, or particles from a simulation box), and project 
them on the sky adding various layers of realism to create a mock survey.

The implemented steps include: 

 1. Volume Remapping (BoxRemap by Carlson and White, see http://mwhite.berkeley.edu/BoxRemap/ )
 2. Translation / Rotation 
 3. Sky Projection (X,Y,Z to RA,DEC,REDSHIFT)
 4. Redshift Distortions (based on peculiar velocities) 
 5. Survey footprint trimming (defined as MANGLE polygon, see http://space.mit.edu/~molly/mangle/ )
 6. Downsample based on sky completeness (defined by region in MANGLE mask) 
 7. Downsample based on radial selection (input file) 

Most of these are optional based on the input configuration, with the sky 
projection being the only mandatory step.

Most functionality is encapsulated into "libraries" that can be used 
independently. The main program `make_survey.cpp` is really just an example 
glue code that uses all libraries.

If you use this code in research that results in publications, please cite the following paper: 

    White, M., Tinker, J., McBride, C.K., 2013, MNRAS, submitted.
    "Mock galaxy catalogs using the quick particle mesh method"
    http://arxiv.org/abs/1309.5532


DEPENDENCIES 
------------

The `make_survey` code has minimal dependencies, requiring only the GNU 
Scientific Library (GSL, see http://www.gnu.org/software/gsl/).

Much of the functionality is included in "libraries", but the source for these 
is included in the codebase (see the `lib/` subdirectory).

The `make_survey` code makes use of the following libraries (bundled):

 - slightly tweaked code from BoxRemap 
		http://mwhite.berkeley.edu/BoxRemap/
 - MANGLE functionality in a minimal (C-only) re-implementation
	  https://github.com/cmcbride/minimal_mangle
 - simple libs and check utilities
	  https://github.com/cmcbride/simple_lib
		https://github.com/cmcbride/check_utils
    

If this functionality is desired, there might be more updated versions at 
those repositories.

The MANGLE functionality is based on a minimal re-implementation (C-only), so 
there is '''no''' dependency on the original MANGLE codebase. The caveat to 
this is that the polygon input is '''not''' flexible.  If the mangle input does 
not conform to the limited versions this code reads, we advise you to use the 
original MANGLE tools to convert to an accepted format.


COMPILATION
-----------

Edit the Makefile to change the compiler and compilation flags. Various 
utilities require both C++ and C compilers. The GSL libraries will be 
automatically detected if the `gsl-config` program is in the default path.

Once the source is downloaded, and GSL is available, one can use GNU make to 
build the code with one of the following:

    % make        # build make_survey main program 
    % make tools  # build additional useful utilities included in source 
    % make all    # what do you think? 


USAGE
-----

    ./make_survey  CONFIG_FILE  MOCK_IN  RDZW_OUT

 `CONFIG_FILE`: Defines most input variables, files, and steps. See `examples/AbacusSummit_c110_params.param` in the source for keyword descriptions that includes the updated w0, wa functionality. For flat lambda CDM, use w0 = -1.0, wa = 0.0 .

 `MOCK_IN`: ASCII input mock file: one object per line with the first 6 columns being positions (3) and velocities (3).
    Requires at least 3 columns of input (redshift space output requires at least 6 columns).
		
 `RDZW_OUT`: ASCII output mock with in 4 columns: `RA DEC REDSHIFT COMP` 
    `COMP` is the sky completeness defined in the input MANGLE mask file (as the weight) 

The general idea is that the config file is written for one general mapping 
(i.e. one set of mocks), and the other two command line arguments can be 
adjusted to easily iterate over many similar realizations.  This allows one to 
easily process a large number of input files using one defining parameter file. 

The `make_survey` code is written to be relatively efficient. It reads the 
input catalog line-by-line (one object at a time), processes each object 
through all required steps, and outputs any object that makes all cuts before 
moving on to the next. This means the whole input catalog is ''never'' fully 
read into memory, nor is the output stored in memory.  Processing speed depends 
on the the required steps (specified in the configuration file) and input 
files. Overall it should be sufficiently fast so that the runtime should not 
require any significant computing resources. One of the longest computational 
costs is the MANGLE mask search, and the code does include a pixelized search 
(only the simple scheme), which can be a dramatic speed-up (and requires a 
pixelized input MANGLE mask). 

For simplicity, ASCII file input and output is currently implemented. If IO 
becomes an issue, the main program (`make_survey.cpp`) can be modified to read 
and write native binary formats. This will both reduce the size of the files 
and speed program access.


CONFIGURATION
-------------

There is a commented parameter file which documents the available keywords:
`examples/dr10_ngc.commented.param` .

Most of these should be relatively easy to understand.

Probably the most confusing configuration options are the remap / translate / 
rotate trio.  There is a non-trivial relation between choosing these values and 
other input such as the simulation volume, the mock survey geometry, mock 
survey depth, and conventions of the sky projection.  There is one additional 
option worth mentioning in this context: `pre-rotate`, but it is independent of 
the others. 

 `pre-rotate`:
    3-element list of integer values of 90 degree rotations around each X, Y, Z-axis of the centered simulation box. 
    (i.e. 2,0,0 = 180 degree rotation around a centered X-axis).  
    This is applied ''before'' remapping, and allows different projections of the box.

 `remap`:
    9-element matrix specifying how the simulation box is remapped into a cuboid.  
    See the http://mwhite.berkeley.edu/BoxRemap/ BoxRemap documentation for full description, 
    and the included `examples/remap.list7.trim.txt` file cached from their genremap utility.

 `translate`:
    3-element translation of X,Y,Z coordinates, applied after remapping but before rotation.

 `rotate`::
    3-element list of angles (in degrees) to rotate around the ''remapped'' and ''translated'' X-, Y-, and Z-axis (in that order).

The sky projection assumes a convention (chosen to match some legacy code), which is defined in 
	`lib/coord_transforms.c`

Picking the best remap / translate / rotate configuration can be a tricky step, and is unlikely to 
have a unique solution.


ADDITIONAL TOOLS
----------------

A number of additional utilities are included in the source. In general, 
running the command without any arguments will output a USAGE statement. 

Brief descriptions of the currently included utilities: 

 `check_survey`:
    Test the `make_survey` input parameter file to see if an INPUT ra/dec/redshift is correctly bounded by the remap / translate / rotate transformation.
    This tool is intended to help choose input parameters for make_survey. 

 `mply_area`:
    Calculate the area in a MANGLE polygon file, optionally trimming by a minimum weight.

 `mply_polyid`:
    Find the POLYID of a MANGLE polygon for an input RA_DEC_FILE.

 `mply_trim`:
    Trim input RA_DEC_FILE based on MANGLE polygon. Can also be run in REVERSE_TRIM mode to process veto masks. 

 `randbox`:
    Create a uniform distribution in X,Y,Z Cartesian coordinates to create randoms (can be fed into make_survey).

 `trim_by_index`:
    Trim input FILE by index assuming the INDEX counts non-comments, non-empty lines in input FILE.
    This can help trim the mock input based on the generated index file with the `make_index 1` configuration option. 

