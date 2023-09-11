import matplotlib.pyplot as plt
import numpy as np
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

def save_blade_data(file_name, data):
    with open(file_name, 'w') as file:
        json.dump(data, file)

def linearize_blade(profiles) -> dict:
    n = len(profiles["0"]["x_cx"])
    z = [profile["z"] for profile in profiles.values()]
    for i in range(n):
        x_cx = []; y_cx = []; x_cv = []; y_cv = []
        for profile in profiles.values():
            x_cx.append(profile.get("x_cx")[i])
            y_cx.append(profile.get("y_cx")[i])
            x_cv.append(profile.get("x_cv")[i])
            y_cv.append(profile.get("y_cv")[i])
        x_cx1 = lstsq(z, x_cx); y_cx1 = lstsq(z, y_cx)
        x_cv1 = lstsq(z, x_cv); y_cv1 = lstsq(z, y_cv)
        for (profile, x1, y1, x2, y2) in zip(profiles.values(), x_cx1, y_cx1, x_cv1, y_cv1):  
            profile["x_cx"][i] = x1
            profile["y_cx"][i] = y1
            profile["x_cv"][i] = x2
            profile["y_cv"][i] = y2
    return profiles

# ------------------------ BEGIN ------------------------ #


#profiles = get_raw_blade_data("../blade_data.json")
#lin_profiles = linearize_blade(profiles)
#save_blade_data("linearized_blade_data.json", lin_profiles)

'''
for profile in profiles.values():
    x = []; y = []
    for point in profile["cx"]:
        x.append(point[0])
        y.append(point[1])
    profile["x_cx"] = x
    profile["y_cx"] = y
    del profile["cx"]
    x = []; y = []
    for point in profile["cv"]:
        x.append(point[0])
        y.append(point[1])
    profile["x_cv"] = x
    profile["y_cv"] = y
    del profile["cx"]
'''

def foo(a: dict) -> dict:
    a[1] = 3
    return a

a = [1,2,3]
b = foo(a)

print(id(a))
print(id(b))

