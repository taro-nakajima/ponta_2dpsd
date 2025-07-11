import sys
import numpy as np
import configparser

args=sys.argv
if (len(args)<2):
	print("usage:  python ponta_2dpsd.py [config file] ")
	sys.exit()

cfg = configparser.ConfigParser()
cfg.read(args[1], encoding='utf-8')

print(f"========== Detector info ===========")
IOfiles = cfg['IO Files']
datalist_file = IOfiles.get('datalist')
output_file = IOfiles.get('outfile')

print(f"========== Detector info ===========")
DetInfo = cfg['Detector Info']
SDD = float(DetInfo.get('SDD'))
A2Center = float(DetInfo.get('A2Center'))
pixelNumX = int(DetInfo.get('pixelNumX'))
pixelNumY = int(DetInfo.get('pixelNumY'))
pixelSizeX = float(DetInfo.get('pixelSizeX'))
pixelSizeY = float(DetInfo.get('pixelSizeY'))
x0 = float(DetInfo.get('centerPixelX'))
y0 = float(DetInfo.get('centerPixelY'))
alpha = float(DetInfo.get('alpha'))
print(f"SDD = {SDD}")
print(f"A2 center = {A2Center}")
print(f"pixelNumX = {pixelNumX}")
print(f"pixelNumY = {pixelNumY}")
print(f"pixelSizeX = {pixelSizeX}")
print(f"pixelSizeY = {pixelSizeY}")
print(f"x0 = {x0}")
print(f"y0 = {y0}")
print(f"alpha = {alpha}")
print("")

print(f"========== Experiment info ===========")
ExpInfo = cfg['Exp Info']
Ei = float(ExpInfo.get('Ei'))
k_len = np.sqrt(Ei/2.072)
print(f"Ei = {Ei} meV (k={k_len} A-1)")
print("")

print(f"========== Experiment info ===========")
SampleInfo = cfg['Sample Info']
C2ofst = float(SampleInfo.get('C2ofst'))
print(f"C2 offset = {C2ofst}")
print("")


print(f"========== Slice info ===========")
SliceInfo = cfg['Slice Info']
Qx_max = float(SliceInfo.get('Qx_max'))
Qx_min = float(SliceInfo.get('Qx_min'))
Qy_max = float(SliceInfo.get('Qy_max'))
Qy_min = float(SliceInfo.get('Qy_min'))
Qz_max = float(SliceInfo.get('Qz_max'))
Qz_min = float(SliceInfo.get('Qz_min'))
axis1 = (SliceInfo.get('Axis1'))
mesh1 = float(SliceInfo.get('Mesh1'))
axis2 = (SliceInfo.get('Axis2'))
mesh2 = float(SliceInfo.get('Mesh2'))
zeroIntFilling = SliceInfo.get('zeroIntFilling')
print(f"Qx_max = {Qx_max}")
print(f"Qx_min = {Qx_min}")
print(f"Qy_max = {Qy_max}")
print(f"Qy_min = {Qy_min}")
print(f"Qz_max = {Qz_max}")
print(f"Qz_min = {Qz_min}")
print(f"Slice plane = {axis1}-{axis2}, ({mesh1} x {mesh2})")
print(f"Zero intensity filling = {zeroIntFilling}")

dQ1=0.0
dQ2=0.0


if (axis1=="x"):
    dQ1=(Qx_max-Qx_min)/mesh1
elif (axis1=="y"):
    dQ1=(Qy_max-Qy_min)/mesh1
elif (axis1=="z"):
    dQ1=(Qz_max-Qz_min)/mesh1

if (axis2=="x"):
    dQ2=(Qx_max-Qx_min)/mesh2
elif (axis2=="y"):
    dQ2=(Qy_max-Qy_min)/mesh2
elif (axis2=="z"):
    dQ2=(Qz_max-Qz_min)/mesh2


## step 0: preparing a matrix with the size of (pixelNumX*pixelNumY,3)

kf_array=np.zeros((pixelNumX*pixelNumY,3))

# step 1: define pixel positions on the yz plane.

for i in range(pixelNumX):
    for j in range(pixelNumY):
        zpos_temp = (float(i)-x0)*pixelSizeX   # Xpixel of the detector -> Z direction
        ypos_temp = (float(j)-y0)*pixelSizeY   # Ypixel of the detector -> Y direction
        pos_vec = np.array([0.0,ypos_temp,zpos_temp])
        pos_vec_norm = pos_vec/np.linalg.norm(pos_vec)
        kf_array[i*pixelNumX+j][0]= pos_vec_norm[0]*k_len       # kf_array[0]=(Xpic=0,Ypic=0), kf_array[1]=(Xpic=0,Ypic=1), kf_array[2]=(Xpic=0,Ypic=2).....
        kf_array[i*pixelNumX+j][1]= pos_vec_norm[1]*k_len
        kf_array[i*pixelNumX+j][2]= pos_vec_norm[2]*k_len

# step 2: z-rotation by alpha

Rot_alpha = np.array( 
    [[ np.cos(np.pi/180.0*(alpha)),  -np.sin(np.pi/180.0*(alpha)),  0 ],
     [ np.sin(np.pi/180.0*(alpha)),  np.cos(np.pi/180.0*(alpha)),  0 ],
     [  0,   0, 1.0 ]])

kf_array = kf_array @ Rot_alpha.T

# step 3: x-translation by SDD

trans_x = np.array([SDD,0,0])

kf_array = kf_array + trans_x

# step 4: z-rotation by A2Center

Rot_A2 = np.array( 
    [[ np.cos(np.pi/180.0*(A2Center)),  -np.sin(np.pi/180.0*(A2Center)),  0 ],
     [ np.sin(np.pi/180.0*(A2Center)),  np.cos(np.pi/180.0*(A2Center)),  0 ],
     [  0,   0, 1.0 ]])

kf_array = kf_array @ Rot_A2.T

# step 5: calculate Q0 vectors (Q-vectors at C2=0) by Q=ki-kf

ki = np.array([k_len,0,0])

Q0=ki-kf_array

Intensity=np.zeros((int(mesh1),int(mesh2)))
SqError=np.zeros((int(mesh1),int(mesh2)))
dataNum=np.zeros((int(mesh1),int(mesh2)))

FH1=open(datalist_file,"r")
for line in FH1:
    temp=line.split()
    intMapFile=temp[0]
    C2=float(temp[1])+C2ofst
    countTime=float(temp[2])
    Rot_C2 = np.array( 
        [[ np.cos(np.pi/180.0*(-C2)),  -np.sin(np.pi/180.0*(-C2)),  0 ],
        [ np.sin(np.pi/180.0*(-C2)),  np.cos(np.pi/180.0*(-C2)),  0 ],
        [  0,   0, 1.0 ]])
    Q_C2rot = Q0 @ Rot_C2.T

    FH2=open(intMapFile,"r")
    for line2 in FH2:
        if "#" not in line2:
            values = line2.split()
            if len(values) == 3:
                Xtemp=int(float(values[0]))
                Ytemp=int(float(values[1]))
                Qx=Q_C2rot[Xtemp*pixelNumX+Ytemp][0]
                Qy=Q_C2rot[Xtemp*pixelNumX+Ytemp][1]
                Qz=Q_C2rot[Xtemp*pixelNumX+Ytemp][2]
                if (Qx_min <= Qx <=Qx_max) and (Qy_min <= Qy <=Qy_max) and (Qz_min <= Qz <=Qz_max):
                    i=0
                    j=0
                    if (axis1=="x"):
                        i = int(((Qx-Qx_min)/dQ1))
                    elif (axis1=="y"):
                        i = int(((Qy-Qy_min)/dQ1))
                    elif (axis1=="z"):
                        i = int(((Qz-Qz_min)/dQ1))

                    if (axis2=="x"):
                        j = int(((Qx-Qx_min)/dQ2))
                    elif (axis2=="y"):
                        j = int(((Qy-Qy_min)/dQ2))
                    elif (axis2=="z"):
                        j = int(((Qz-Qz_min)/dQ2))

                    Intensity[i][j]+=float(values[2])/countTime
                    SqError[i][j]+=float(values[2])/countTime**2.0
                    dataNum[i][j]+=1


FHR=open(output_file,"w")
FHR.write(f"#Q{axis1}  Q{axis2}  Intensity  Error  dataNum\n")
for i in range(int(mesh1)):
    for j in range(int(mesh2)):
        Q1=0.0
        Q2=0.0
        if (axis1 == "x"):
            Q1=Qx_min+dQ1*float(i)
        elif (axis1 == "y"):
            Q1=Qy_min+dQ1*float(i)
        elif (axis1 == "z"):
            Q1=Qz_min+dQ1*float(i)

        if (axis2 == "x"):
            Q2=Qx_min+dQ2*float(j)
        elif (axis2 == "y"):
            Q2=Qy_min+dQ2*float(j)
        elif (axis2 == "z"):
            Q2=Qz_min+dQ2*float(j)

        if Intensity[i][j] > 0:
            FHR.write("{0}  {1}  {2}  {3}  {4}\n".format(Q1,Q2,Intensity[i][j],np.sqrt(SqError[i][j]),dataNum[i][j]))
        else:
            FHR.write("{0}  {1}  {2}  {3}  {4}\n".format(Q1,Q2,zeroIntFilling,0,dataNum[i][j]))

    FHR.write("\n")

FHR.close()
