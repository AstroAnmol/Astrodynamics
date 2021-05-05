# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# %%
data=pd.read_csv("DirectMars_0_file.csv")
##
Earth=pd.read_csv("Results/Earth_orbit_file.csv")
Mars=pd.read_csv("Results/Mars_orbit_file.csv")
# %%
E_R1=np.array(Earth['Radius_1 (km)'])
E_R2=np.array(Earth['Radius_2 (km)'])
E_R3=np.array(Earth['Radius_3 (km)'])
# %%
M_R1=np.array(Mars['Radius_1 (km)'])
M_R2=np.array(Mars['Radius_2 (km)'])
M_R3=np.array(Mars['Radius_3 (km)'])
# %%
R1=np.array(data['Radius_1'][:-4])
R2=np.array(data['Radius_2'][:-4])
R3=np.array(data['Radius_3'][:-4])
# Velocity
V1=np.array(data['Velocity_1'][:-4])
V2=np.array(data['Velocity_2'][:-4])
V3=np.array(data['Velocity_3'][:-4])
#mass
m=np.array(data['Mass'][:-4])
#Adjoint 
AR1=np.array(data['ARadius_1'][:-4])
AR2=np.array(data['ARadius_2'][:-4])
AR3=np.array(data['ARadius_3'][:-4])
AV1=np.array(data['AVelocity_1'][:-4])
AV2=np.array(data['AVelocity_2'][:-4])
AV3=np.array(data['AVelocity_3'][:-4])
Am=np.array(data['AMass'][:-4])
# Time
time=np.array(data['Time'][:-4])
# Thrust
Thrust=np.array(data['Thrust'][:-4])
# %%
# Calculate the SF
AV=np.sqrt(AV1*AV1 + AV2*AV2 + AV3*AV3)
# Switching
SF=AV/m - Am
SF=np.array(SF)
R1_thrusting=[]
R2_thrusting=[]
R3_thrusting=[]
for i in range(0,SF.size):
    if SF[i]>0:
        R1_thrusting.append(R1[i])
        R2_thrusting.append(R2[i])
        R3_thrusting.append(R3[i])
    elif SF[i]<0:
        print(i)
R1_thrusting=np.array(R1_thrusting)
R2_thrusting=np.array(R2_thrusting)
R3_thrusting=np.array(R3_thrusting)
print(SF)
# %%
fig =plt.figure(1)
plt.plot(E_R1,E_R2,label='Earth')
plt.plot(M_R1,M_R2,label='Mars')
plt.plot(R1,R2, ls='-.', label='Spacecraft Coasting')
plt.plot(R1_thrusting[:656],R2_thrusting[:656],label='Spacecraft Thrusting 01')
plt.plot(R1_thrusting[657:],R2_thrusting[657:],label='Spacecraft Thrusting 02')
plt.legend()
plt.axis('square')
plt.xlabel('x-HCI')
plt.ylabel('y-HCI')
plt.show()
# %%
fig =plt.figure(2)
ax = fig.gca(projection='3d')
ax.set_box_aspect((np.ptp(R1), np.ptp(R2), np.ptp(R3)))
ax.plot(E_R1,E_R2,E_R3,label='Earth')
ax.plot(M_R1,M_R2,M_R3,label='Mars')
ax.plot(R1,R2,R3,label='Spacecraft')
plt.legend()
plt.show()
# %%
plt.figure(3)
plt.plot(range(0,SF.size), SF)
plt.show()