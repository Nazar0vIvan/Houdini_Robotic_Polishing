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
    return x, k*x+b

profiles = get_raw_blade_data("../blade_data.json")

# spl = scipy.interpolate.CubicSpline(x, y)
# x = np.linspace(min(x), max(x), num = 50)
# plot(x, spl(x))

z = []
for profile in profiles.values():
    z.append(profile["z"])

n = len(profiles["0"]["convex"])
for i in range(n):
    y = []
    for profile in profiles.values():
        y.append(profile.get("convex")[i][1])
    x, y_new = lstsq(z, y)
    for (profile, y) in zip(profiles.values(), y_new):
        profile["convex"][i][1] = y

    # plt.plot(xs, ys, 'o-', linewidth = 1)
    # plt.plot(z, y, 'o-', linewidth = 1)

# with open("blade_ruled_surface.json", 'w') as file:
#     json.dump(profiles, file)

for i in range(n):
    y = []
    for profile in profiles.values():
        y.append(profile.get("convex")[i][1])
    plt.plot(z, y, '-', linewidth=1)
    plt.grid(True)

plt.show()












