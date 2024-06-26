# Numerical Computation of seismograms in 3D plane layered medium
### Last version May 2024 fixes an amplitude normalization problem introduced few versions ago
### Results are now checked againt theoretical response, see EXAMPLE/check_amplitude.ipynb 

Axitra Force & Moment versions 

Should you encounter any problem using the program, please feel free to contact me at:
Olivier Coutant: olivier.coutant@univ-grenoble-alpes.fr

Thank you to cite the following paper if you publish results obtained with axitra:
F. Cotton and Coutant O., 1997, Dynamic stress variations due to shear faults in a plane-layered medium, 
GEOPHYSICAL JOURNAL INTERNATIONAL,Vol 128, 676-688

DOI: 10.5281/zenodo.4603722

TO COMPILE
==========

In the Makefile 
1) select the appropriate options for the compiler (FC) and its options (FFLAGS)
2) the code uses openmp parallelization
3) then
make clean  #clean the object, modules
make axitra #compile the greeen's function program
make convms #compile the convolution program
make python #compile python3.7 wrappers

To run
4) edit the "<header>.data" and "<header>.hist" files according to your needs
   <header> is either 'axi' or 'axi_suffix' where suffix is an optional parameters passed as positional argument to the two programs
5) run the codes:
./axitra [suffix]
./convms [suffix]
6) the seismograms are written as sac files

To run under python, see "example.py" or jupyter notebook "example.ipynb" examples under EXAMPLE directory


NEW FEATURES
============

This "new" version includes the following changes
- a python interface and wrappers to the code
- all input/output files can be given a "suffix" in order to have several configurations in the same directory
  axitra and convms are then run with the suffix as argument
- the namelist is reduced to a few parameters, nr,ns and nc are dynamically counted
- The duration (tl) and periodicity (xl) can be automatically computed if needed when the values are set to (<=0)
in the axi.data namelist

and previous changes

- porting to fortran95 allows for dynamic allocations
- openMP directives for use on shared memory multicore computer
- force sources and moment sources versions are now available in two separate programs (and directories)
 with almost the same functionnalities: multiple source; multiple receiver; same convolution process
- possibility to enter position in lat,lon coordinates
  when latlon variable is set to .true., the program will output coordinates
  in cartesian coordinates in .xyz files. The origin is arbibrary, these cartesian coordinates
  can then be used as input to the program after setting latlon variable to .false.

- possibility to enter an "axi.date" file that gives an origin date for synthetic
  when wrinting in sac format

- attenuation formulation follows Martin Galis iand Emmanuel Chaljub suggestion. 2 mechanisms
Futterman and Kjartansson.


USING AXITRA
============

The reflectivity code consist of 2 executables: 
(1) axitra computes the Green's function in the frequency domain for 
- 5 double couples + isotropic source, at each receiver and for each source (MOMENT VERSION)
- 3 unidirectionanal forces (FORCE VERSION)

(2) convms convolves these Green's function
by the appropriate source time function's, and for the chosen mechanism. The 
output are the seismograms in the time domain in sac format.

In this last version, preprocessing has been removed. axitra now by default
- run in double precision
- uses geographical XYZ axes (No radial, tranverse,...)
- uses harwell routine to compute bessel functions


axitra is then ready to run. The different data files are:
(1) input for axitra:
- "axi[.suffix].data" +  file giving station coordinates + file giving sources 
  coordinates
(2) output from axitra:
- axi[.suffix].head (information about convergence)
- axi[.suffix].res  (binary transfert function)
(3) input for the convolution program convms
- axi[.suffix].hist

PARAMETERS (file "parameter.f90")
==========
- nkmax= If you have chose to store the bessel functions (which do not depend 
  upon the frequency), nkmax values will be stored for the first nkmax radial 
  wavenumber, if the number of iteration is larger than nkmax, the bessel 
  functions are recomputed each time for those nkmax+... iterations. This 
  value depends on the available memory and the convergence rate. 
  You can try to run the null frequency only to have a rough idea of
  how long it takes. Set nkmax to 3/2*number_of_iteration_for_freq=0 should be 
  ok if you don't compute too high frequencies (the higher the frequency, the 
  slower the convergence). (ps: You really save a lot of time by storing the 
  bessel functions)
- nkmin= minimum number fo iteration for the wavenumber summation

The parameters defined in "axi.data" are:
(1) in the namelist:
- nfreq=number of frequency to compute
- tl=time length of the seismogram (sec), if set to <=0, tl is computed by axitra to contain all arrivals
- aw=coefficient for frequency imaginary part, omega=(2*pi*freq, aw*pi/tl)
- xl=medium periodicity (m or km, be consistent with the other variables), if set to <=0, automatically computed
- ikmax=max number of iteration, generate an error message if overpassed
- latlon=coordinates given in (lat,lon) (.true.) or in m/km (.false.)
  when latlon is true, distance unit is set to meter (depth). velocities are in m/s and density in kg/m^3
- freesurface= free surface is set at Z=0 (.true.) or upper space is infinite (.false.)
- sourcefile=source description file
- statfile=receiver description file

(2) in free format and for each layer:
- thickness (or depth of the upper interface), vp, vs, rho, qp,qs
  Units are [m, m/s and Kg/m^3], or [km, km/s and Kg/km^3]
  If you specify the upper interface depth, it must be set to 0. for the
  free surface. !!!! rho unit must be consistent with the length unit (m or km) !!!!

INPUT DESCRIPTION FILES + AXI.HIST
==================================
in free format, 

I) For the receivers:
receiver_index,x,y,z
 or
receiver_index,lat,lon,z

II) For the sources:
source_index, x, y, z
 or
source_index, lat, lon, z 

III) For the file AXI.HIST

MOMENT VERSION:
 2 possibilities depending on the data you chose to
give: moment or displacement+width+length:
source_index, moment, strike, dip, rake, 0., 0., t0
  or
source_index, disp , strike, dip, rake, width, length, t0

  If you want a null contribution for a given source, just set 
  the corresponding  moment (resp disp) to zero

You can compute seismograms generated by an explosive source by setting a 
negative surface equal to -1 (either negative width or length). The explosion moment
is then read as the given displacement ('disp' variable)

FORCE VERSION:
source_index, fx_amplitude, fy_amplitude, fz_amplitude, total_amplitude, time_delay

AXIS CONVENTION
===============
The convention for the axis is:
x=north, y=east, z=upward.
!!!!!!!! Warling !!!!!!! This convention is different in the computation where the vertical axis is 
positive downward. 

When given in geographical coordinates, axitra read x=latitude then y=longitude

LIMITATIONS/NOTES
=================
-Since the programs uses an horizontal wavenumber decomposition,
the only geometry not allowed is to put the receivers and the sources 
at the same depth. 
-The convergence rate depends in part on the depth difference between source 
  and receiver, the closer they are, the longer the convergence is.
- In some (rare) cases, the program may generate an underflow error when 
  computing the term "exp(-i*gamma*(z-zref)) where z may be
  the receiver or source depth, zref  the source or layer depth. This
  may occurs when you compute very high frequencies (KHz). A solution
 in that case is to incorporate a new artificial layer between z and zref
to reduce the (z-zref) term.
- An artificial vertical plane wave may be (rarely) seen when the source is
close to the free surface and when you look at regional distances. This
wave is seen only on the terms involving the J0 bessel function. It 
decreases very quickly as a function of the source depth but has no
geometrical spreading (it is a plane wave) and may be visible at the
beginning of the seismogram. There is nothing yet to do to remove
it completly. Using a recipy given by Bob Herrmann, the horizontal
wavenumbers are given a small fractional part, this reduces greatly
the effet of this spurious plane wave.

At this time and after several years, we do not know of any injuries due to 
the use of this program, anyway use it at your own risk...

