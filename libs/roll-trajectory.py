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

# profiles = get_raw_blade_data("../blade_data.json")

for profile in profiles.values():
    x = []; y = []
    for point in profile["convex"]:
        x.append(point[0])
        y.append(point[1])
    profile["x_cx"] = x
    profile["y_cx"] = y
    del profile["convex"]
    x = []; y = []
    for point in profile["concave"]:
        x.append(point[0])
        y.append(point[1])
    profile["x_cv"] = x
    profile["y_cv"] = y
    del profile["concave"]

save_blade_data("new_blade_data_newformat.json")

# spl = scipy.interpolate.CubicSpline(x, y)
# x = np.linspace(min(x), max(x), num = 50)
# plot(x, spl(x))












