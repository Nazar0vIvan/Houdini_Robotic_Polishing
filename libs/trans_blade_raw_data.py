
import json

def parse_json(filePath: str) -> dict:
    try:
        file = open(filePath, 'r')
        try:
            return trans_blade_raw_data(json.load(file))
        except json.decoder.JSONDecodeError as error:
            print("Decoding JSON has failed: " + str(error))  
    except OSError as error:
        print(error)

def trans_blade_raw_data(airfoil: dict) -> dict:
    for profile in airfoil.values():
        z = profile['z']

        convex_points = [] 
        for i,(xc, yc) in enumerate(zip(profile['xc'], profile['yc'])):
            convex_points.append([xc,yc])
        profile["convex"] = convex_points

        concave_points = []
        for i,(xk, yk) in enumerate(zip(profile['xk'], profile['yk'])):
            concave_points.append([xk,yk])
        profile["concave"] = concave_points
            
        del profile['xk']; del profile['yk']; del profile['xc']; del profile['yc']
    return airfoil

print(parse_json("242s_raw.json"))

