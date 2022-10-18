# Written by Jameel Ahmed Syed
# Email id: j.syed@innopolis.university

# Simple transformations

from sympy import Matrix, diff, det, solve, atan2, simplify, trigsimp, sqrt
import numpy as np
from traceback import format_exc
import datetime

def setup_symbolic(symbol):
    global symbolic
    symbolic = symbol

def rotz(t):
    """Rotates the frame about Z Axis"""
    if symbolic is True:
        from sympy import cos, sin
    else:
        from numpy import cos, sin
    rot = Matrix([[cos(t), -sin(t), 0, 0],
                  [sin(t), cos(t),  0, 0],
                  [0,       0,      1, 0],
                  [0,       0,      0, 1]])
    return rot

def rotx(t):
    """Rotates the frame about X Axis"""
    if symbolic is True:
        from sympy import cos, sin
    else:
        from numpy import cos, sin
    rot = Matrix([[1,       0,      0, 0],
                  [0, cos(t), -sin(t), 0],
                  [0, sin(t), cos(t),  0],
                  [0,   0,      0,     1]])
    return rot


def roty(t):
    """Rotates the frame about Y Axis"""
    if symbolic is True:
        from sympy import cos, sin
    else:
        from numpy import cos, sin
    rot = Matrix([[cos(t),  0, sin(t), 0],
                  [0,       1,   0,    0],
                  [-sin(t), 0, cos(t), 0],
                  [0,       0,   0,    1]])
    return rot


def tranx(d):
    """Translates the frame on X Axis"""
    transx = Matrix([[1, 0, 0, d],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
    return transx

def trany(d):
    """Translates the frame on Y Axis"""
    transy = Matrix([[1, 0, 0, 0],
                     [0, 1, 0, d],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
    return transy

def tranz(d):
    """Translates the frame on Z Axis"""
    transz = Matrix([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, d],
                     [0, 0, 0, 1]])
    return transz

def print_matrix(input_matrix, name_matrix='T'):
    """Prints a Matrix in a clear Matrix format
    Just a little fancy stuff to make the matrix printing readable"""
    try:
        a = np.array(input_matrix)
        if a.ndim > 0:
            print(name_matrix, '=')
            s = [[str(e) for e in row] for row in a]
            lens = [max(map(len, col)) for col in zip(*s)]
            fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
            table = [fmt.format(*row) for row in s]
            print('\n'.join(table), '\n')
        else:
            print(name_matrix, '= ', end="")
            print(a)

    except Exception:
        error_file = open("error_logs.txt", "a")
        error_file.writelines(f"{datetime.datetime.now()}\n{format_exc()}\n\n")
        error_file.close()