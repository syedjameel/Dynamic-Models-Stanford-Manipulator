# Written by Jameel Ahmed Syed
# Email id: j.syed@innopolis.university

# Forward, Inverse, Jacobian for the Stanford Manipulator with spherical wrist

from SimpleTranformations import *
symbol_names = ['t1', 't2', 't3', 't4', 't5', 't6', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6']


def forward(joint_angles, link_lengths, euler_wrist, symbolic):
    """Forward Kinematics Solution For Stanford Manipulator
    with Spherical Wrist in any valid euler_wrist arrangement
    :param joint_angles: Contains the joint angles t1, t2, t4, t5, t6 and d3(d3 is not an angle its displacement)"""
    if symbolic is False:
        t1, t2, t3, t4, t5, t6 = joint_angles
        d1, d2, d3, d4, d5, d6 = link_lengths

    if symbolic is True:
        from sympy.solvers import solve
        from sympy import symbols
        t1 = symbols('t1')
        t2 = symbols('t2')
        t4 = symbols('t4')
        t5 = symbols('t5')
        t6 = symbols('t6')
        d2 = symbols('d2')
        d3 = symbols('d3')
        d6 = symbols('d6')

    setup_symbolic(symbolic)
    from SimpleTranformations import rotz, rotx, tranx, tranz
    global t01, t12, t23, t34, t45, t56, t03, t36, t06, p0w, p06, X, Y, Z, t02, t03, t04, t05
    t01 = rotz(t1)                                  # Rotation around Z axis
    t12 = rotx(t2) * tranx(d2)                      # Rotation around Z axis and Translation along X axis
    t23 = tranz(d3)                                 # Translation along Z axis
    t34 = get_rot(euler_wrist[0])(t4)               # Rotation around Z axis    #rotz(t4)
    t45 = get_rot(euler_wrist[1])(t5)               # Rotation around X axis    #rotx(t5)
    t56 = get_rot(euler_wrist[2])(t6) * tranz(d6)   # Rotation around Z axis and Translation along Z axis   #rotz(t6)
    t03 = t01 * t12 * t23                           # Transformation from 0 to 3
    t36 = t34 * t45 * t56                           # Transformation from 3 to 6
    t06 = t03 * t36                                 # Transformation from 0 to 6 (Complete Forward Kinematic Solution)
    # t06md = t01 * t12 * t23 * t34 * t45 * rotz(t6)  # Transformation from 0 to Wrist by eliminating the d6 translation
    p06 = t06[0:3, 3]                       # This is the P0 to wrist position using Pipers Method
    p0w = p06 - d6 * t06[0:3, 2]
    t02 = t01 * t12
    t03 = t02 * t23
    t04 = t03 * t34
    t05 = t04 * t45
    X = t06[0, 3]
    Y = t06[1, 3]
    Z = t06[2, 3]
    # if symbolic is True:
    #    d = diff(t06[2, 3], d3)
    #    print('Diff = ', d)
    T = {'T01': t01, 'T12': t12, 'T23': t23, 'T34': t34, 'T45': t45, 'T56': t56, 'T03': t03, 'T36' : t36, 'T06': t06, 'P0w' : p0w, 'P06': p06, 'X' : X, 'Y': Y, 'Z' : Z, 'T02': t02, 'T04': t04, 'T05': t05}
    return T

# ToDo
def inverse(position, wrist_configuration):
    pass

def jacobian_geometrical(tranforamtions):
    """Returns the Jacobian Matrix from the Forward kinematics Geometrically"""
    j = ['j1', 'j2', 'j3', 'j4', 'j5', 'j6']

    for var in j:
        globals()[var] = Matrix([[0], [0], [0], [0], [0], [0]])

    z0 = Matrix([[0], [0], [1]])
    p0 = Matrix([[0], [0], [0]])
    j1 = z0.cross((p06-p0)).row_insert(3, z0)
    j2 = t01[0:3, 0].cross((p06-t01[0:3, 3])).row_insert(3, t01[0:3, 0])
    j3 = t02[0:3, 2].row_insert(3, Matrix([[0], [0], [0]]))
    j4 = t03[0:3, 2].cross((p06-t03[0:3, 3])).row_insert(3, t03[0:3, 2])
    j5 = t04[0:3, 0].cross((p06-t04[0:3, 3])).row_insert(3, t04[0:3, 0])
    j6 = t05[0:3, 2].cross((p06-t05[0:3, 3])).row_insert(3, t05[0:3, 2]) # if the position is the same as p0w then change this thing to t04[0:3, 3]

    cols = [j1, j2, j3, j4, j5, j6]

    J = Matrix([])
    for j in range(len(cols)):
        J = J.col_insert(j, cols[j])
    #J = {'J1': j1, 'J2': j2, 'J3': j3, 'J4': j4, 'J5': j5, 'J6': j6}

    return J

def jacobian_numerically(function_list, t1, t2, d3, t4, t5, t6, d2, d6):
    """Returns the Jacobian Matrix Numerically
    :param function_list: It contains the equations of Position and Orientation
            i.e. X, Y, Z, ThetaX, ThetaY, ThetaZ
    :param t1: Angle around the Joint 1
    :param t2: Angle around the Joint 2
    :param d3: Distance along the Joint 3
    :param t4: Angle around the Joint 4
    :param t5: Angle around the Joint 5
    :param t6: Angle around the Joint 6"""
    from math import cos, sin
    q_original = []
    t = [t1, t2, d3, t4, t5, t6]
    for i in t:
        q_original.append(i)
    d = q_original.copy()
    j = []
    J = Matrix([])
    for fun in range(len(function_list)):
        for i in range(len(t)):
            fx = eval(str(function_list[fun]))
            dh = 1.0e-8
            t[i] = d[i] + dh
            t1, t2, d3, t4, t5, t6 = t
            fxh = eval(str(function_list[fun]))
            f = (fxh - fx) / dh
            t1, t2, d3, t4, t5, t6 = update(d)
            t = update(d)
            j.append(f)
        J.col_insert(fun-1, Matrix([j]))
    return j

def update(d):
    a = []
    a = d
    return a
def get_rot(axis):
    """Gets the rotations as per the Euler's wrist configurations"""
    if axis == 'x' or axis == 'X':
        rot = rotx
    elif axis == 'y' or axis == 'Y':
        rot = roty
    elif axis == 'z' or axis == 'Z':
        rot = rotz
    return rot

# ToDo
def check_valid_wrist_configuration():
    pass