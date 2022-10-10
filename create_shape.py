#from ase import Atoms
from ase.io import read, write
from ase.build import make_supercell
import numpy as np
import pandas as pd
from scipy.integrate import quad
#from scipy.spatial.transform import Rotation
#from scipy import special
#from ase.optimize import BFGS
#from ase.calculators.emt import EMT
import argparse
import sys
import math
import matplotlib.pyplot as plt
from scipy.interpolate import splrep,splev
# from scipy.signal import argrelextrema


parser = argparse.ArgumentParser(description="Create a structure of a 2D layer which follows a parametric function")
parser.add_argument("--func", type=int, default=2, help="The parametric function will be (by default) \
                    1 - sine, 2 - ellipsoidic, 3 - nanotube, 4 - hypotrochoid 5- two sines 6- numeric ")
parser.add_argument("--dir", type=int, default=3, help="Create the wrinkle along the specified lattice vector (default=the smaller of 1 or 2. vector)")
parser.add_argument("--A", type=float, required=False, default=0, help="Amplitude of the wave or the 1st radius")
parser.add_argument("--length", type=float, required=False, default=0,help="Wavelength of the wave or the 2nd radius")
parser.add_argument("--C", type=float, default=1.0, help="For some curves a 3rd parameter is needed.")
parser.add_argument("--struct", type=str, required=True, help="The unit cell of the 2D layer, properly centered along z!")
parser.add_argument("--intype", type=str, default='xyz', help="Type of input file as used in ASE (default=xyz)")
parser.add_argument("--adjust", type=int, default=1, help="Adjust the wavelength, 0, or the amplitude, 1 (default=1)")
parser.add_argument("--out", type=str, help="Output file (default=L'length'A'A'.in)")
parser.add_argument("--outtype", type=str, default='aims', help="Type of output file as used in ASE (default=aims)")
parser.add_argument("--unit", type=int, default=0, required=False, help="for wrinkles, this option forces the number of unit cells")
parser.add_argument("--chiral", type=int, nargs=2, default=0, required=False, help="chiral structure rolling vector, for the moment tested for hexagonal cell")


args = parser.parse_args()

#
# m=5
# n=1
outfile='geometry_supercell.in'
#

def chiral_unitcell(inputfile,outfile,m,n):

    ################# reading the geometry
    single_cell = read(inputfile, format='aims')
    gcd_value = math.gcd(m,n)
    # n=n/gcd_value
    # m=m/gcd_value
    # formula for hexagonal unit cell
    l1=-(2*n+m)
    l2=2*m+n
    P=[[m,n,0],[l1,l2,0],[0,0,1]]
    final_cell = make_supercell(single_cell,P,wrap=False)  
    final_cell.rotate(np.degrees(np.arctan(-final_cell.get_cell()[0][1]/final_cell.get_cell()[0][0])),(0,0,1),rotate_cell=True)
    cm=final_cell.get_center_of_mass()
    final_cell.translate([0,0,-cm[2]])
    # positions = final_cell.get_positions()
    # positions[:,2] = 0 # make things on the ground
    # final_cell.set_positions(positions)
    write(outfile,final_cell, format='aims',wrap=False)
    return final_cell,l1,l2

    
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
elif (args.func == 5):
  def pfunction(t, a, b, c): # two sines function with seven parameters (including one shift) a0*sin(a1*x+a2)+a3*sin(a4*s+a5)+a6
    a0=a*(-1.08971)
    a1=1 #I had to change the a1 and a4 to get the perfect periodicity, other values are one of the fitted data
    a2=-3.17079
    a3=a*0.110569
    a4=3
    a5=-6.23166
    a6=0
    x = t/2/np.pi*b
    y = a0*np.sin(a1*t+a2)+a3*np.sin(a4*t+a5)+a6
    dx = b/2/np.pi
    dy = a0*a1*np.cos(a1*t+a2)+a3*a4*np.cos(a4*t+a5)
    return(x, y, dx, dy)
  def arclength(x1, x2, a, b, c):
    a0=a*(-1.08971)
    a1=0.987038
    a2=-3.17079
    a3=a*0.110569
    a4=2.91407
    a5=-6.23166
    return quad(lambda t: np.sqrt(math.pow(a0*a1*np.cos(a1*t+a2)+a3*a4*np.cos(a4*t+a5),2)+math.pow(b/2/np.pi,2)), x1, x2)[0]
elif (args.func == 6):
    # file = open('wrinkle_profile.txt') 
    filename='wrinkle_profile.txt'
    data = pd.read_csv(filename, sep="\s+", header=None)
    # lines = file.readlines()
    positions=data.iloc[3:124,2:4] #selecting the wrinkle+MoS2 
    positions = positions.replace({',':'.'},regex=True).apply(pd.to_numeric,1)
    positions_A =positions*10**10 #change of unit to Angstrom
    x=np.array(positions_A.iloc[::,0])
    y=np.array(positions_A.iloc[::,1])
    # x_m= list(map( lambda t: ((2*np.pi)/(max(positions_A.iloc[::,0])-min(positions_A.iloc[::,0])))*t, positions_A.iloc[::,0]))
    # y_m= list(map( lambda t: (2/(max(positions_A.iloc[::,1])-min(positions_A.iloc[::,1])))*t, positions_A.iloc[::,1]))
    Bspline=  splrep(positions_A.iloc[::,0],positions_A.iloc[::,1] ,k=4)
    New_K=400
    t =np.linspace(0,9500,New_K)
    plt.plot(t,splev(t,Bspline),'o')
    # plt.xlim([0,8200])
    dy=[]
    for i in range(len(t)-1): dy.append((splev(t,Bspline)[i+1]-splev(t,Bspline)[i])/(t[i+1]-t[i]))
    ds=[]
    for i in range(len(dy)): ds.append(np.sqrt(1+dy[i]**2)*(t[2]-t[1]))

    arclength_curve= sum(ds)
    #
    curvature=[]
    dyy=[]
    for i in range(len(t)-2):
      # dyy.append((dy[i+1]-dy[i])/((t[i+1]-t[i])))
        dyy.append((splev(t,Bspline)[i+2]-2*splev(t,Bspline)[i+1]+splev(t,Bspline)[i])/((t[i+1]-t[i])**2))

    for i in range(len(t)-2):
        curvature.append(dyy[i]/((1+dy[i]**2)**(3/2)))

else:
  print(" [fatal error] parametric function ",args.func, " not defined")
  sys.exit()

wrinkle = args.dir
chiral = args.chiral
if chiral!=0: m,n =args.chiral
###message section, for making the code neat I put all the messages about the inputs here
if (wrinkle > 1 or wrinkle <0 and chiral!=0):
    print( "[WARNING] wrinkle direction ", args.dir, " not defined") 
    print("[WARNING] wrinkle direction is set to None")
    # wrinkle =None
unit = args.unit
if unit == 0:
    print( '[WARNING] the number of unit cells in a the wrinkle is not present therefore it is calculated from length and amplitud') 
if unit != 0 and args.func == 3:
    print( '[fatal error] I cannot always keep the number of unit cells and make a prefect nanotube') #later I will correct for nanotube as I have to change in both direction, and comment later for the user
    sys.exit()   
if (args.A == 0 or args.length == 0):
        print('[WARNING] no amplitude or wavelength is not given if it is a unitcell only case there will be no problem')
if (args.A == 0 or args.length == 0) and (args.func != 2 and args.func !=6 and args.chiral ==0):
        sys.exit('[fatal error ] for the selected shapes I need wavelength and amplitude')
if (args.A == 0 and args.length == 0 and args.func == 2 and args.unit == 0 and args.chiral ==0):   
        sys.exit('[fatal error] for making a cycloidic structure without amplitude and wavelength I need the number of unit cells if the structure is chiral provide rolling vector')
####################################
inputfile = args.struct

if chiral !=0:
    final_cell,l1,l2 =chiral_unitcell(inputfile, outfile, m, n)
    single_cell = final_cell
else:
     single_cell=read(inputfile, format=args.intype)
#single_cell.wrap(pbc=[1,1,1],center=(0.0,0.0,0.0))
# print(single_cell.get_positions())

vec1=single_cell.get_cell()[0]
vec2=single_cell.get_cell()[1]
vec3=single_cell.get_cell()[2]

a=single_cell.get_cell_lengths_and_angles()[0]
b=single_cell.get_cell_lengths_and_angles()[1]

if chiral!=0:
    wrinkle=0
elif wrinkle==None:
  if a > b:
    wrinkle=1
  else:
    wrinkle=0

if wrinkle==0:
  print("[info] Wrinkle along lattice direction 1!", a)
  nowrinkle=1
else:
  print("[info] Wrinkle along lattice direction 2!", b)
  nowrinkle=0

alat=single_cell.get_cell_lengths_and_angles()[wrinkle]
func = args.func
#Amplitude and wavelength of the wrinkle in Angstrom
if (args.A ==0 and args.length == 0 and chiral ==0):
    amplitude=alat*args.unit/(2*np.pi)
    wavelength=4*amplitude
elif(chiral ==0 and func == 2): 
     amplitude=alat*args.unit/2/(2*np.pi)
elif (chiral !=0 and args.A==0 and func == 2):
    wavelength = np.linalg.norm(vec1)
    amplitude = 0.25* wavelength
elif (chiral !=0 and func == 3):
    wavelength = np.linalg.norm(vec1)/2/(2*np.pi)
    amplitude = wavelength
elif (chiral !=0 and args.A !=0):
    wavelength = np.linalg.norm(vec1)
    amplitude = args.A
else:
   amplitude=args.A
   wavelength=args.length


cparam=args.C


#Now, should wavelength (0) or amplitude (1) be adjusted?
adjusting=args.adjust
#Note that for the hypotrochoid, C will always be adepted

if (args.func == 4):
  multipi = np.lcm(np.floor(amplitude).astype(int),np.floor(wavelength).astype(int))/amplitude
  print("multipi ",multipi)
else:
  multipi = 1

#circumference=wavelength*special.ellipe(np.sqrt(1-math.pow(amplitude/(wavelength/4),2)))
if (args.func == 6):
    circumference= arclength_curve
elif (chiral!=0):
    circumference=np.linalg.norm(vec1)
else:
    circumference=arclength(0,2*np.pi*multipi,amplitude,wavelength,cparam)

print("[info] Initial wavelength ", wavelength)
print("[info] Initial amplitude ", amplitude)
print("[info] Arclength of the function ", circumference)
if unit == 0: 
    howoften=round(circumference / alat)
else:
        howoften = unit
print("[info] How often does the cell fits into this (rounded) ",howoften)

length=howoften*alat
print("[info] The resulting length of the path would be ",length)

if (args.func == 6):
    if (length ==arclength_curve):
        print("final length is determined")
    else:
        # x_m= list(map( lambda t: ((2*np.pi)/(max(positions_A.iloc[::,0])-min(positions_A.iloc[::,0])))*t, positions_A.iloc[::,0]))
        # y_m= list(map( lambda t: (2/(max(positions_A.iloc[::,1])-min(positions_A.iloc[::,1])))*t, positions_A.iloc[::,1]))
        length

if abs(length-circumference)==0:
  print("[info] initial is final") 
else:
  if length > circumference:
    print("[info] we need to increase")
    factor=1.0
  else:
    print("[info] we need to reduce")
    factor=-1.0
  wavetest=wavelength
  amptest=amplitude
 # print('amptest after amptest', amptest)
  ctest=cparam
  circumtest=circumference
  if (unit == 0 and args.func !=6):
      while factor*(circumtest-length) < 0:
         
        if (args.func == 3):
          wavetest = wavetest + factor*0.00001*wavelength
          amptest = amptest + factor*0.00001*amplitude
        elif (args.func == 4):
          ctest = ctest + factor*0.0001*cparam
        elif (args.func == 2):
          amptest= amptest + factor*0.0001*amplitude
          #print('amptest in func 2',amptest)
          #print('factor*(circumtest-length) in func 2', factor*(circumtest-length))
        else:
          if adjusting==0:
            wavetest = wavetest + factor*0.00001*wavelength
          if adjusting==1:
            amptest = amptest + factor*0.00001*amplitude
        circumtest = arclength(0,2*np.pi*multipi,amptest,wavetest,ctest)
        if wavetest < 0 or amptest < 0 or ctest < 0:
          print("wavetest, amptest ",wavetest,amptest)
          sys.exit('Wavelength or amplitude negative! Reasonable wave form cannot be found.')
  if (args.func ==6):
      a=1
  else: 
  #    print('amptest in the else ',amptest)
       while factor*(circumtest-length) < 0:
          
       #  print('amptest in the while loop ',amptest)
       #  print('while status',factor*(circumtest-length),'amptest',amptest)
         amptest = amptest + factor*0.00001*amplitude 
         circumtest = arclength(0,2*np.pi*multipi,amptest,wavetest,ctest)
         if wavetest < 0 or amptest < 0 or ctest < 0:
          print("amptest ",amptest)
          sys.exit('Wavelength or amplitude negative! Reasonable wave form cannot be found.')
  howoften=round(circumtest / alat)
  length=howoften*alat
  wavelength=wavetest
  amplitude=amptest
  cparam=ctest
  circumference=circumtest
  print("[info] Arclength of the function", circumference)
  #print(arclength(0,2*np.pi*multipi,amplitude,wavelength,cparam))
  print('[info] number of unit cell is',howoften)


final_cell = single_cell.copy()
del final_cell[:]
if (10*amplitude <100 or 10*amplitude > 300):
       print('[WARNING] check your cell dimension for vacuum size!')
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

if wrinkle == 1:
  rotation_axis=np.cross(vec3,vec1)
#  rotation_axis=vec2
elif chiral!=0:
    rotation_axis=vec2
else:
  rotation_axis=np.cross(vec3,vec2)
#  rotation_axis=-vec1

rotation_axis=rotation_axis/np.linalg.norm(rotation_axis)
print('[info] rotation axis:',rotation_axis)

stepping = 0.001 * np.pi
t0 = 0
i = 0
counting=1

while counting <= howoften:
  pathlength=arclength(t0,t0+i*stepping,amplitude,wavelength,cparam)
#  print("initial pathlength",pathlength)
  tmp_cell = None
  tmp_cell = single_cell.copy()
  for natom in range(len(tmp_cell.positions)):
   # print("initial pos",tmp_cell.positions[natom])
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
