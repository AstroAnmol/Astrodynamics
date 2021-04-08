# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# %%
data=pd.read_csv("Trial__file.csv")

# %%
R1=np.array(data['Radius_1'][:-4])
R2=np.array(data[' Radius_2'][:-4])
R3=np.array(data[' Radius_3'][:-4])
# Velocity
V1=np.array(data[' Velocity_1'][:-4])
V2=np.array(data[' Velocity_2'][:-4])
V3=np.array(data[' Velocity_3'][:-4])
#mass
m=np.array(data[' Mass'][:-4])
#Adjoint 
AR1=np.array(data[' ARadius_1'][:-4])
AR2=np.array(data[' ARadius_2'][:-4])
AR3=np.array(data[' ARadius_3'][:-4])
AV1=np.array(data[' AVelocity_1'][:-4])
AV2=np.array(data[' AVelocity_2'][:-4])
AV3=np.array(data[' AVelocity_3'][:-4])
Am=np.array(data[' AMass'][:-4])
# Time
time=np.array([data[' Time'][:-4]])
# %%
# Calculate the SF
AV=np.sqrt(AV1*AV1 + AV2*AV2 + AV3*AV3)
# Switching
SF=AV/m - Am
SF=np.array(SF)
R1_thrusting=[]
R2_thrusting=[]
R3_thrusting=[]
Thrust=np.zeros(SF.size)
for i in range(0,SF.size):
    if SF[i]>0:
        Thrust[i]=1
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
ax = fig.gca(projection='3d')
ax.plot(R1,R2,R3)
plt.show()
# %%
fig =plt.figure(2)
#plt.plot(R1,R2, ls='-.')
plt.plot(R1_thrusting,R2_thrusting)
plt.show()
# %%
plt.plot(m, SF)
plt.show()