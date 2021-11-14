import platform
import tkinter
import numpy as np
#from numpy.polynomial.polynomial import Polynomial, polyder, polyfromroots, polyval
#import kp
#import pyshader as ps
#from scipy import optimize
import PIL
from PIL import ImageDraw, ImageTk
import math
import time
import multiprocessing as mp
from multiprocessing import Pool
from tkinter import *

COLORS = [[140,95,102],[172,188,165],[232,185,171],[224,152,145],[203,118,158],[44, 17, 30],[49, 33, 35],[78, 94, 69]]  # Fall scheme
#COLORS = [[157,34,18],[94,201,160],[90,137,218],[20,44,86],[184,158,90],[38,8,4],[36,72,73]] # Red/Blue scheme
NO_ROOT = [0,0,0]
# 4k = 3840 x 2160
MAX_X = 800
MAX_Y = 600
MAX_THREADS = 8  #mp.cpu_count()
SCALE = 1
OFFSET_X = 0
OFFSET_Y = 0
TOLERANCE = .001
MAX_ITER = 200

iterations = [[0] * MAX_Y for i in range(MAX_X) ]

class Poly:
    def __init__(self, terms):
        self.terms = []
        for t in terms:
            self.terms.append([t[0], t[1]])

    def eval(self, point):
        sum = complex(0,0)
        for t in self.terms:
            sum += (point ** t[1]) * t[0]
        return sum

    def newton(self, deriv, point, iter):
        numer = self.eval(point)
        denom = deriv.eval(point)
        if denom == 0:
            raise ZeroDivisionError
        newPoint = point - (numer / denom)
        if abs(newPoint - point) < TOLERANCE:
            root = newPoint
        else:
            (root, iter) = self.newton(deriv, newPoint, iter + 1)
        return root, iter

    def toString(self):
        string = ""
        for t in self.terms:
            if t[0] != 0:
                string += " + " + str(normalize(t[0])) + "x^" + str(t[1])
        return string

def normalize(point):
    newPoint = complex(round(point.real, 3), round(point.imag, 3))
    return newPoint

def fillSerial():
    for x in range(MAX_X):
        for y in range(MAX_Y):
            # Generate point to calcuate, centered at (ORIGIN_X,ORIGIN_Y), scaled in by SCALE
            point = complex((-(MAX_X/2)+x+OFFSET_X)*SCALE, (-(MAX_Y/2)+y+OFFSET_Y)*SCALE)
            iter = 0
            try:
                (root, iter) = poly.newton(deriv, point, iter)
                root = normalize(root)
                grid[y][x] = COLORS[roots.index(root)]
            except ZeroDivisionError:
                grid[y][x] = NO_ROOT
            for rgbVal in range(3):
                grid[y][x][rgbVal] = max(0, min(255, grid[y][x][rgbVal] - (3 * iter)))

def fillParallel(id, p):    
    startY = int((MAX_Y * id) / p)
    endY = int((MAX_Y * (id + 1)) / p)
    gridPart = []
    for x in range(MAX_X):
        for y in range(startY, endY):
            # Generate point to calcuate, centered at (ORIGIN_X,ORIGIN_Y), scaled in by SCALE
            point = complex((-(MAX_X/2)+x+OFFSET_X)*SCALE, (-(MAX_Y/2)+y+OFFSET_Y)*SCALE)
            iter = 0
            try:
                (root, iter) = poly.newton(deriv, point, iter)
                root = normalize(root)
                #iterations[x][y] = iter
                gridPart.append([COLORS[roots.index(root)], iter])
            except ZeroDivisionError:
                gridPart.append([NO_ROOT, iter])
    return gridPart

roots = []
# Function 1: x^3 - 1
roots.append(normalize(complex(1, 0)))
roots.append(normalize(complex(-.5, math.sqrt(3)/2)))
roots.append(normalize(complex(-.5, -math.sqrt(3)/2)))

# Function 2: x^8 + 15x^4 -16
#roots.append(normalize(complex(1,0)))
#roots.append(normalize(complex(-1,0)))
#roots.append(normalize(complex(0,1)))
#roots.append(normalize(complex(0,-1)))
#roots.append(normalize(complex(math.sqrt(2),math.sqrt(2))))
#roots.append(normalize(complex(-math.sqrt(2),-math.sqrt(2))))
#roots.append(normalize(complex(-math.sqrt(2),math.sqrt(2))))
#roots.append(normalize(complex(math.sqrt(2),-math.sqrt(2))))

poly = Poly([[1,3],[-1,0]])
deriv = Poly([[3,2]])

grid = np.ndarray([MAX_Y,MAX_X,3], dtype=np.uint8)

drawType = "Parallel"
numT = 1

startTime = time.perf_counter()

if drawType == "Serial":
    fillSerial()
elif drawType == "Parallel":
    #almost working
    numT = MAX_THREADS
    results = []
    if __name__ == '__main__':
        with Pool(processes=numT) as pool:
            for i in range(numT):
                res = pool.apply_async(fillParallel, (i, numT,))
                results.insert(i, res)
            for i in range(numT):
                val = results[i].get()
                index = 0
                for x in range(MAX_X):
                    startY = int((MAX_Y * i) / numT)
                    endY = int((MAX_Y * (i + 1)) / numT)
                    for y in range(startY, endY):
                        grid[y][x] = val[index][0]
                        for rgbVal in range(3):
                            grid[y][x][rgbVal] = max(0, min(255, grid[y][x][rgbVal] - (3 * val[index][1])))
                        index += 1
            pool.close()
            pool.join()
#elif drawType == "GPU":
    # Creates Kompute Manager
    #mgr = kp.Manager()

    # Create Kompute Tensors to hold data
    #tensor_in_a = kp.Tensor([2,2,2])
    #tensor_in_b = kp.Tensor([1,2,3])
    #tensor_out = kp.Tensor(grid)

if __name__ == "__main__":
    endTime = time.perf_counter()

    bitmap = PIL.Image.fromarray(grid)
    #bitmap = bitmap.filter(ImageFilter.DETAIL)

    draw = ImageDraw.Draw(bitmap)
    
    textLabel = drawType + " " + str(MAX_X) + "x" + str(MAX_Y) + " in Python"
    textLabel += "\n" + platform.system() + " " + platform.machine() + " Using " + str(numT) + "/" + str(mp.cpu_count()) + " Threads"
    #if drawType == "GPU":
    #    textLabel += mgr.get_device_properties()["device_name"]
    textLabel += "\nScale: " + str(SCALE) + "  Offset: " + str(OFFSET_X) + "," + str(OFFSET_Y)
    textLabel += "\nFunction: " + poly.toString()
    textLabel += "\nIn " + str(round((endTime-startTime)*1000,2)) + " ms"
    draw.text((2,2), textLabel)

    #bitmap.show()

    root = Tk()
    root.geometry(str(MAX_X) + "x" + str(MAX_Y))
    frame = Frame(root)
    frame.pack()

    myImage = ImageTk.PhotoImage(bitmap)

    label = tkinter.Label(image=myImage)
    label.image = myImage
    label.pack()

    root.title("Newton's Fractal")
    root.mainloop()


"""
# Initialize Kompute Tensors in GPU
mgr.eval_tensor_create_def([tensor_in_a,tensor_in_b,tensor_out])

@ps.python2shader
def compute_shader_multiply(index=("input", "GlobalInvocationId", ps.ivec3),
                            data1=("buffer", 0, ps.Array(ps.f32)),
                            data2=("buffer", 1, ps.Array(ps.f32)),
                            data3=("buffer", 2, ps.Array(ps.f32))):
    i = index.x # Fetch the current run index being processed
    data3[i] = data1[i] * data2[i]

mgr.eval_algo_data_def(
    [tensor_in_a, tensor_in_b, tensor_out],
    compute_shader_multiply.to_spirv())

mgr.eval_tensor_sync_local_def([tensor_out])

print(tensor_out.data())
"""
