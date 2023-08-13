import matplotlib.pyplot as plt
import numpy as np
import scipy
import json

def get_raw_blade_data(filePath: str) -> dict:
    try:
        file = open(filePath, 'r')
        try:
            return json.load(file)
        except json.decoder.JSONDecodeError as error:
            print("Decoding JSON has failed: " + str(error))  
    except OSError as error:
        print(error)

def plot(x, y, xlabel = "x", ylabel = "y"): 
    plt.plot(x, y, '-', linewidth=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()

def lstsq(x, y):
    x = np.array(x)
    y = np.array(y)
    A = np.vstack([x, np.ones(len(x))]).T
    k, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return k*x+b

def save_blade_data(file_name):
    with open(file_name, 'w') as file:
        json.dump(profiles, file)
    file.close()

profiles = get_raw_blade_data("../blade_data.json")

z = []
for profile in profiles.values():
    z.append(profile["z"])

n = len(profiles["0"]["convex"])

for i in range(n):
    x_cx = []; y_cx = []; x_cv = []; y_cv = []
    for profile in profiles.values():
        x_cx.append(profile.get("convex")[i][0])
        y_cx.append(profile.get("convex")[i][1])
        x_cv.append(profile.get("concave")[i][0])
        y_cv.append(profile.get("concave")[i][1])
    x_cx1 = lstsq(z, x_cx)
    y_cx1 = lstsq(z, y_cx)
    x_cv1 = lstsq(z, x_cv)
    y_cv1 = lstsq(z, y_cv)
    for (profile, x1, y1, x2, y2) in zip(profiles.values(), x_cx1, y_cx1, x_cv1, y_cv1):
        profile["convex"][i][0] = x1
        profile["convex"][i][1] = y1
        profile["concave"][i][0] = x2
        profile["concave"][i][1] = y2

# save_blade_data("new_blade_data.json")

# spl = scipy.interpolate.CubicSpline(x, y)
# x = np.linspace(min(x), max(x), num = 50)
# plot(x, spl(x))


""" a = -2.05915523e-08 
b = -7.53749861e-02 
c = 6.34916362e+00

x = np.arange(10)

plot(x, a*x**2 + b*x + c) """

""" 
for i in range(n):
    x = []
    for profile in profiles.values():
        x.append(profile.get("convex")[i][1])
    plt.plot(z, x, '-', linewidth=1)
    plt.grid(True)

plt.show() 
"""












