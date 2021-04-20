#from ase import Atoms
from ase.io import read, write
import numpy as np
from scipy.integrate import quad
#from scipy.spatial.transform import Rotation
#from scipy import special
#from ase.optimize import BFGS
#from ase.calculators.emt import EMT
import argparse
import sys
import math

parser = argparse.ArgumentParser(description="Create a structure of a 2D layer which follows a parametric function")
parser.add_argument("--func", type=int, default=1, help="The parametric function will be (by default) 1 - sine, 2 - ellipsoidic, 3 - nanotube, 4 - hypotrochoid")
parser.add_argument("--dir", type=int, help="Create the wrinkle along the specified lattice vector (default=the smaller of 1. or 2. vector)")
parser.add_argument("--A", type=float, required=True, help="Amplitude of the wave or the 1st radius")
parser.add_argument("--length", type=float, required=True, help="Wavelength of the wave or the 2nd radius")
parser.add_argument("--C", type=float, default=1.0, help="For some curves a 3rd parameter is needed.")
parser.add_argument("--struct", type=str, required=True, help="The unit cell of the 2D layer, properly centered along z!")
parser.add_argument("--intype", type=str, default='xyz', help="Type of input file as used in ASE (default=xyz)")
parser.add_argument("--adjust", type=int, default=1, help="Adjust the wavelength, 0, or the amplitude, 1 (default=1)")
parser.add_argument("--out", type=str, help="Output file (default=L'length'A'A'.in)")
parser.add_argument("--outtype", type=str, default='aims', help="Type of output file as used in ASE (default=aims)")

args = parser.parse_args()

#

if (args.func == 1):
   def pfunction(t, a, b, c): # the sine-like parametric function. x=t/2pi * wavelength, y=amplitude * sin(t)
    x = t/2/np.pi*b
    y = a*np.sin(t)
    dx = b/2/np.pi
    dy = a*np.cos(t)
    return (x, y, dx, dy)
   def arclength(x1, x2, a, b, c):
    return quad(lambda t: np.sqrt(math.pow(a*np.cos(t),2)+math.pow(b/2/np.pi,2)), x1, x2)[0]
   print('the sine-like parametric function. x=t/2pi * wavelength, y=amplitude * sin(t) is requested!')  
elif (args.func == 2):    # the elliptic-like parametric function.
  def pfunction(t, a, b, c): # first the upper half, shifted such that x goes from 0 till wavelength/2
    if (t <= np.pi):
      x = 1/4 * b - b/4 * np.cos(t)
      y = a * np.sin(t)
      dx = b/4 * np.sin(t)
      dy = a * np.cos(t)
    else:                 # the second (lower) half
      x = 3/4 * b + b/4 * np.cos(t)
      y = a * np.sin(t)
      dx = -b/4 * np.sin(t) # might be that we have to change this in order to have a continuous
      dy = a * np.cos(t)    # and sensible normal vector, n=(dy, -dx), which always points to the
    return(x, y, dx, dy)    # "upper" side of the 2D material
  def arclength(x1, x2, a, b, c):
    return quad(lambda t: np.sqrt(math.pow(a*np.cos(t),2)+math.pow(b/4*np.sin(t),2)), x1, x2)[0]
  print('the elliptic-like parametric function wrinkle is requested!')
elif (args.func == 3):
  def pfunction(t, a, b, c): # the circle-like parametric function. x=cos(t) * wavelength, y=amplitude * sin(t)
    x = b*np.cos(t)
    y = a*np.sin(t)
    dx = -b*np.sin(t)
    dy = a*np.cos(t)
    return (x, y, dx, dy)
  def arclength(x1, x2, a, b, c):
    return quad(lambda t: np.sqrt(math.pow(a*np.cos(t),2)+math.pow(b*np.sin(t),2)), x1, x2)[0]
  print('the circle-like parametric function (i.e. nanotube). x=cos(t) * wavelength, y=amplitude * sin(t) is requested')
elif (args.func == 4):
  def pfunction(t, a, b, c): # Hypotrochoid function, a=R, b=r, c=d
    x = (a - b)*np.cos(t) + c*np.cos((a - b)/b*t)
    y = (a - b)*np.sin(t) - c*np.sin((a - b)/b*t)
    dx = ((b - a)*(c*np.sin(t*(a/b - 1)) + b*np.sin(t)))/b
    dy = ((b - a)*(c*np.cos(t*(a/b - 1)) - b*np.cos(t)))/b
    return(x, y, dx, dy)
  def arclength(x1, x2, a, b, c):
    return quad(lambda t: np.sqrt(math.pow(((b - a)*(c*np.sin(t*(a/b - 1)) + b*np.sin(t)))/b,2)+math.pow(((b - a)*(c*np.cos(t*(a/b - 1)) - b*np.cos(t)))/b,2)), x1, x2)[0]
else:
  print("parametric function ",args.func, " not defined")
  sys.exit()

wrinkle = args.dir
if wrinkle > 1 or wrinkle <0:
    print( "[WARNING] wrinkle direction ", args.dir, " not defined") 
    print("[WARNING] wrinkle direction is set to None")
    wrinkle =None
    
inputfile = args.struct

single_cell=read(inputfile, format=args.intype)
single_cell.wrap(pbc=[1,1,0])
#single_cell.wrap(pbc=[1,1,1],center=(0.0,0.0,0.0))
print(single_cell.get_positions())

vec1=single_cell.get_cell()[0]
vec2=single_cell.get_cell()[1]
vec3=single_cell.get_cell()[2]

a=single_cell.get_cell_lengths_and_angles()[0]
b=single_cell.get_cell_lengths_and_angles()[1]


if wrinkle==None:
  if a > b:
    wrinkle=1
  else:
    wrinkle=0

if wrinkle==0:
  print("Wrinkle along lattice direction 1!", a)
  nowrinkle=1
else:
  print("Wrinkle along lattice direction 2!", b)
  nowrinkle=0

#Amplitude and wavelength of the wrinkle in Angstrom
amplitude=args.A
wavelength=args.length
cparam=args.C

alat=single_cell.get_cell_lengths_and_angles()[wrinkle]

#Now, should wavelength (0) or amplitude (1) be adjusted?
adjusting=args.adjust
#Note that for the hypotrochoid, C will always be adepted

if (args.func == 4):
  multipi = np.lcm(np.floor(amplitude).astype(int),np.floor(wavelength).astype(int))/amplitude
  print("multipi ",multipi)
else:
  multipi = 1

#circumference=wavelength*special.ellipe(np.sqrt(1-math.pow(amplitude/(wavelength/4),2)))
circumference=arclength(0,2*np.pi*multipi,amplitude,wavelength,cparam)

print("Initial wavelength ", wavelength)
print("Initial amplitude ", amplitude)
print("Arclength of the function ", circumference)
howoften=round(circumference / alat)
print("How often does the cell fits into this (rounded) ",howoften)
length=howoften*alat
print("The resulting length of the path would be ",length)

if abs(length-circumference)==0:
  print("initial is final")
else:
  if length > circumference:
    print("we need to increase")
    factor=1.0
  else:
    print("we need to reduce")
    factor=-1.0
  wavetest=wavelength
  amptest=amplitude
  ctest=cparam
  circumtest=circumference
  while factor*(circumtest-length) < 0:
     
    if (args.func == 3):
      wavetest = wavetest + factor*0.00001*wavelength
      amptest = amptest + factor*0.00001*amplitude
    elif (args.func == 4):
      ctest = ctest + factor*0.0001*cparam
    elif (args.func == 2):
      amptest= amptest + factor*0.0001*amplitude
     #print(amptest)
    else:
      if adjusting==0:
        wavetest = wavetest + factor*0.00001*wavelength
      if adjusting==1:
        amptest = amptest + factor*0.00001*amplitude
    circumtest = arclength(0,2*np.pi*multipi,amptest,wavetest,ctest)
    if wavetest < 0 or amptest < 0 or ctest < 0:
      print("wavetest, amptest ",wavetest,amptest)
      sys.exit('Wavelength or amplitude negative! Reasonable wave form cannot be found.')
  print("final wavelength and amplitude: ",wavetest,amptest)
  howoften=round(circumtest / alat)
  length=howoften*alat
  wavelength=wavetest
  amplitude=amptest
  cparam=ctest
  circumference=circumtest
  print("Arclength of the function", circumference)
  print(arclength(0,2*np.pi*multipi,amplitude,wavelength,cparam))
  print(howoften)


final_cell = single_cell.copy()
del final_cell[:]
if (10*amplitude <100 or 10*amplitude > 300):
       print('[WARNING] check your cell dimension!')
if (args.func == 3 or args.func == 4):
  if wrinkle == 0:
    final_cell.set_cell([(10*amplitude,0,0),vec2,(0,0,10*amplitude)])
  else:
    final_cell.set_cell([vec1,(0,10*amplitude,0),(0,0,10*amplitude)])
else:
  if wrinkle == 0:
    final_cell.set_cell([(wavelength/a)*vec1,vec2,(0,0,10*amplitude)])
  else:
    final_cell.set_cell([vec1,(wavelength/b)*vec2,(0,0,10*amplitude)])

###
### Note for myself: it is not correct to rotate around one of the cell vectors!!
### we have to rotate around the vector perpendicular to the plane defined by the
### wrinkle vector and the z axis!!! Rotation around cell vectors only works if
### the corresponding wrinkle vector and the z axis build a plane perpendicular
### to the other vector, i.e., for rectangular grids...
###

if wrinkle == 0:
  rotation_axis=np.cross(vec3,vec1)
#  rotation_axis=vec2
else:
  rotation_axis=np.cross(vec3,vec2)
#  rotation_axis=-vec1

rotation_axis=rotation_axis/np.linalg.norm(rotation_axis)
print(rotation_axis)

stepping = 0.00001 * np.pi
t0 = 0
i = 0
counting=1

while counting <= howoften:
  pathlength=arclength(t0,t0+i*stepping,amplitude,wavelength,cparam)
#  print("initial pathlength",pathlength)
  tmp_cell = None
  tmp_cell = single_cell.copy()
  for natom in range(len(tmp_cell.positions)):
    print("initial pos",tmp_cell.positions[natom])
    while tmp_cell.positions[natom][wrinkle] > pathlength:
      i += 1
      pathlength=arclength(t0,t0+i*stepping,amplitude,wavelength,cparam)
    x, y, dx, dy = pfunction(t0+i*stepping,amplitude,wavelength,cparam)
#    print("x,y,dx,dy",x,y,dx,dy)
    n=[dy,-dx]
    n=n/np.linalg.norm(n)*tmp_cell.positions[natom][2]
#    print("n",n)
    tmp_cell.positions[natom][wrinkle] = x + n[0]
    tmp_cell.positions[natom][2] = y + n[1]
#    print("final pos",tmp_cell.positions[natom])
    i=0
    pathlength=arclength(t0,t0+i*stepping,amplitude,wavelength,cparam)
#    print("pathlength",pathlength)
  final_cell.extend(tmp_cell)
  while pathlength < alat:
    i += 1
    pathlength=arclength(t0,t0+i*stepping,amplitude,wavelength,cparam)
  t0=t0+i*stepping
  i=0
  print("t0 ",(t0/np.pi)," pi")
  counting += 1

#print(final_cell.positions)

if args.out==None:
  outfile='L' + str(wavelength) + 'A' + str(amplitude) + '.in'
else:
  outfile=args.out

#write('test.xyz',final_cell, vec_cell=True)
write(outfile,final_cell, format=args.outtype)
