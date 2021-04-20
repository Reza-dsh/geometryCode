
the general tool to make structures using parametric functions 
implemented functions are wrinkles in form of ellipsoid and sine wave as well as nanotubes with parametric function.


usage: create_shape.py [-h] [--func FUNC] [--dir DIR] --A A --length LENGTH
                       [--C C] --struct STRUCT [--intype INTYPE]
                       [--adjust ADJUST] [--out OUT] [--outtype OUTTYPE]

Create a structure of a 2D layer which follows a parametric function

optional arguments:
  -h, --help         show this help message and exit
  --func FUNC        The parametric function will be (by default) 1 - sine, 2
                     - ellipsoidic, 3 - nanotube, 4 - hypotrochoid
  --dir DIR          Create the wrinkle along the specified lattice vector
                     (default=the smaller of 1. or 2. vector)
  --A A              Amplitude of the wave or the 1st radius
  --length LENGTH    Wavelength of the wave or the 2nd radius
  --C C              For some curves a 3rd parameter is needed.
  --struct STRUCT    The unit cell of the 2D layer, properly centered along z!
  --intype INTYPE    Type of input file as used in ASE (default=xyz)
  --adjust ADJUST    Adjust the wavelength, 0, or the amplitude, 1 (default=1)
  --out OUT          Output file (default=L'length'A'A'.in)
  --outtype OUTTYPE  Type of output file as used in ASE (default=aims)



