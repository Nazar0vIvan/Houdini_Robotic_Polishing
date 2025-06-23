import numpy as np
import json
import matplotlib.pyplot as plt

# with open("99.01.25.242.json", 'r') as file:
#     blade = json.load(file)

# n = len(blade)
# x = []; y = []; z = []
# for i in range(n-2):
#     x.append(blade[i]["cx"][5][0])
#     y.append(blade[i]["cx"][5][1])
#     z.append(blade[i]["cx"][5][2])


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.set_xlabel('$X$', fontsize=20)
# ax.set_ylabel('$Y$', fontsize=20)
# ax.set_zlabel('$Z$', fontsize=20)

# # Data for a three-dimensional line
# ax.plot3D(x, y, z, 'bo-')
# plt.show()

with open("242_short.json", 'r') as file:
    blade = json.load(file)

x = []; y = []; z = []
for pf in blade.values():
    x.append(pf["x_cx"][4])
    y.append(pf["y_cx"][4])
    z.append(pf["z"])
    print(pf["x_cx"][4])
    print(pf["y_cx"][4])
    print(pf["z"])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel('$Z$', fontsize=20)
ax.set_ylabel('$X$', fontsize=20)
ax.set_zlabel('$Y$', fontsize=20)

# Data for a three-dimensional line
ax.plot3D(z, x, y, 'bo-')
plt.show()