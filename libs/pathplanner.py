import numpy as np
from math import floor

class PathPlanner: 
    def __init__(self, fps: int, path: list = None):
        self.__fps = fps
        self.__path = [] if path == None else path

    @property
    def path(self):
        return self.__path
    
    def add_lin(self, start: list, stop: list, feed: float, is_const_orien = False) -> None:
        dist = np.linalg.norm(np.array(stop[0:3]) - np.array(start[0:3]))
        extra_points_count = floor(dist/feed*self.__fps)
        if(is_const_orien):
            lin_path = np.linspace(start[0:3], stop[0:3], extra_points_count).tolist()
            lin_path.pop()
            lin_path = [p.extend(start[3:6]) for p in lin_path]
        else:
            lin_path = np.linspace(start, stop, extra_points_count).tolist()
            lin_path.pop()

        self.__path.extend(lin_path)

    def add_path(self, path: list) -> None:
        self.__path.extend(path)

        
def save_path(path: list):
    with open('path.txt', 'w') as file:
        file.writelines("[" + ", ".join("{:6.3f}".format(el) for el in pos) + "]\n" for pos in path)
    file.close()

def path2cls(feed, rpm, path = []):
    with open("tool_trajectory.cls", 'w') as file:
        file.write(f"PAINT/SPEED,{feed*60/1000}\n")
        file.write(f"SPINDL/RPM,{rpm},CLW\n")
        for p in path:
            str_p = ','.join(f'{el:.4f}' for el in p)
            file.write(f"GOTO/{str_p}\n")