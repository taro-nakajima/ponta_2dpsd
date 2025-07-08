import numpy as np
import configparser

cfg = configparser.ConfigParser()
cfg.read('config.ini', encoding='utf-8')

DetInfo = cfg['Detector Info']
SDD = float(DetInfo.get('SDD'))
A2Center = float(DetInfo.get('A2Center'))
pixelNumX = int(DetInfo.get('pixelNumX'))
pixelNumY = int(DetInfo.get('pixelNumY'))
pixelSizeX = float(DetInfo.get('pixelSizeX'))
pixelSizeY = float(DetInfo.get('pixelSizeY'))
x0 = float(DetInfo.get('centerPixelX'))
y0 = float(DetInfo.get('centerPixelY'))

print(f"SDD = {SDD}")
print(f"A2 center = {A2Center}")
print(f"pixelNumX = {pixelNumX}")
print(f"pixelNumY = {pixelNumY}")
print(f"pixelSizeX = {pixelSizeX}")
print(f"pixelSizeY = {pixelSizeY}")
print(f"x0 = {x0}")
print(f"y0 = {y0}")

ExpInfo = cfg['Exp Info']
wavelen = float(ExpInfo.get('wavelen'))
print(f"wavelen = {wavelen}")


SampleInfo = cfg['Sample Info']
C2ofst = float(SampleInfo.get('C2ofst'))
print(f"C2 offset = {C2ofst}")

SliceInfo = cfg['Slice Info']
Qx_max = float(SliceInfo.get('Qx_max'))
Qx_min = float(SliceInfo.get('Qx_min'))
Qy_max = float(SliceInfo.get('Qy_max'))
Qy_min = float(SliceInfo.get('Qy_min'))
Qz_max = float(SliceInfo.get('Qz_max'))
Qz_min = float(SliceInfo.get('Qz_min'))
mesh = float(SliceInfo.get('mesh'))
zeroIntFilling = SliceInfo.get('zeroIntFilling')

print(f"Qx_max = {Qx_max}")
print(f"Qx_min = {Qx_min}")
print(f"Qy_max = {Qy_max}")
print(f"Qy_min = {Qy_min}")
print(f"Qz_max = {Qz_max}")
print(f"Qz_min = {Qz_min}")

dQx=(Qx_max-Qx_min)/mesh
dQy=(Qy_max-Qy_min)/mesh
dQz=(Qz_max-Qz_min)/mesh


k_len = 2.0*np.pi/wavelen

kf_array=np.zeros((pixelNumX*pixelNumY,3))

R = np.array( 
    [[ np.cos(np.pi/180.0*(A2Center)),  -np.sin(np.pi/180.0*(A2Center)),  0 ],
     [ np.sin(np.pi/180.0*(A2Center)),  np.cos(np.pi/180.0*(A2Center)),  0 ],
     [  0,   0, 1.0 ]])

for i in range(pixelNumX):
    for j in range(pixelNumY):
        zpos_temp = (float(i)-x0)*pixelSizeX   # Xpixel of the detector -> Z direction
        ypos_temp = (float(j)-y0)*pixelSizeY   # Ypixel of the detector -> Y direction
        pos_vec = np.array([SDD,ypos_temp,zpos_temp])
        pos_vec_norm = pos_vec/np.linalg.norm(pos_vec)
        kf_array[i*pixelNumX+j][0]= pos_vec_norm[0]*k_len       # kf_array[0]=(Xpic=0,Ypic=0), kf_array[1]=(Xpic=0,Ypic=1), kf_array[2]=(Xpic=0,Ypic=2).....
        kf_array[i*pixelNumX+j][1]= pos_vec_norm[1]*k_len
        kf_array[i*pixelNumX+j][2]= pos_vec_norm[2]*k_len



kf_array = kf_array @ R.T

ki = np.array([k_len,0,0])

Q0=ki-kf_array

Intensity=np.zeros((int(mesh),int(mesh)))
SqError=np.zeros((int(mesh),int(mesh)))
dataNum=np.zeros((int(mesh),int(mesh)))


FH1=open("datalist.txt","r")
for line in FH1:
    temp=line.split()
    intMapFile=temp[0]
    C2=float(temp[1])+C2ofst
    countTime=float(temp[2])
    C2rot = np.array( 
        [[ np.cos(np.pi/180.0*(-C2)),  -np.sin(np.pi/180.0*(-C2)),  0 ],
        [ np.sin(np.pi/180.0*(-C2)),  np.cos(np.pi/180.0*(-C2)),  0 ],
        [  0,   0, 1.0 ]])
    Q_C2rot = Q0 @ C2rot.T

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
                    i_Qx = int(((Qx-Qx_min)/dQx))
                    j_Qy = int(((Qy-Qy_min)/dQy))
                    Intensity[i_Qx][j_Qy]+=float(values[2])/countTime
                    SqError[i_Qx][j_Qy]+=float(values[2])/countTime**2.0
                    dataNum[i_Qx][j_Qy]+=1


FHR=open("SlicedData.txt","w")

for i in range(int(mesh)):
    for j in range(int(mesh)):
        Qx=Qx_min+dQx*i
        Qy=Qy_min+dQy*j
        if Intensity[i][j] > 0:
            FHR.write("{0}  {1}  {2}  {3}  {4}\n".format(Qx,Qy,Intensity[i][j],np.sqrt(SqError[i][j]),dataNum[i][j]))
        else:
            FHR.write("{0}  {1}  {2}  {3}  {4}\n".format(Qx,Qy,zeroIntFilling,0,dataNum[i][j]))

    FHR.write("\n")

FHR.close()
