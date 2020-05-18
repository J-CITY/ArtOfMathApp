import numpy as np
from PIL import Image
from enum import Enum

class Format(Enum):
    RGB = 1
    HEX = 2
    BGR = 3
    
TYPE = Format.BGR
IMAGE = "1.png"
PALETTRE_SIZE = 5

def clamp(x): 
  return max(0, min(x, 255))

def palette(img):
    arr = np.asarray(img)
    palette, index = np.unique(asvoid(arr).ravel(), return_inverse=True)
    palette = palette.view(arr.dtype).reshape(-1, arr.shape[-1])
    count = np.bincount(index)
    order = np.argsort(count)
    return palette[order[::-1]]

def asvoid(arr):
    arr = np.ascontiguousarray(arr)
    return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[-1])))

def printRGB(palettre):
    print('[', end='')
    for i, p in enumerate(palettre):
        print(p, end='')
        if (i != len(palettre)-1):
            print(',', end='')
    print(']')

def printBGR(palettre):
    print('[', end='')
    for i, p in enumerate(palettre):
        print(p[::-1], end='')
        if (i != len(palettre)-1):
            print(',', end='')
    print(']')


def printHEX(palettre):
    print('[', end='')
    for i, p in enumerate(palettre):
        print("#{0:02x}{1:02x}{2:02x}".format(clamp(p[0]), clamp(p[1]), clamp(p[2])), end='')
        if (i != len(palettre)-1):
            print(',', end='')
    print(']')

img = Image.open(IMAGE, 'r').convert('RGB')
p = palette(img)

if (TYPE == Format.RGB) :
    printRGB(p[:PALETTRE_SIZE])
elif (TYPE == Format.BGR) :
    printBGR(p[:PALETTRE_SIZE])
else:
    printHEX(p[:PALETTRE_SIZE])

