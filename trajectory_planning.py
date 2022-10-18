import numpy as np
from art import tprint
from matplotlib.pyplot import *
import time
from numpy import pi
from sympy import Matrix
import numpy as np
from traceback import format_exc
import datetime
from numpy import cos, sin

j = []                  # [j1, j2, j3, j4, j5, .... so on till NUM_JOINTS]
j_copy = []
NUM_JOINTS = 6          # Number of joints
dt = 0.05               # 50 ms Controller frequency

j1 = [0, 90, 5, 500/253]      # Joint 1 = [q0, qf, velocity_max, acceleration_max]
j2 = [0, 60, 5, 500/253]      # Joint 2 = [q0, qf, velocity_max, acceleration_max]
j3 = [0, 30, 5, 500/253]      # Joint 3 = [q0, qf, velocity_max, acceleration_max]
j4 = [0, 40, 5, 500/253]      # Joint 4 = [q0, qf, velocity_max, acceleration_max]
j5 = [0, 30, 5, 500/253]      # Joint 5 = [q0, qf, velocity_max, acceleration_max]
j6 = [0, 30, 5, 500/253]      # Joint 6 = [q0, qf, velocity_max, acceleration_max]

# Automatically adding all the j1, j2, j3 ... lists to j till NUM_JOINTS
num_joints = list()
for d in range(NUM_JOINTS):
    num_joints.append(str(int(d + 1)))
    te = tuple(num_joints)
for k, v in list(globals().items()):
    if k.startswith('j') and k.endswith(te):
        j.append(v.copy())
        j_copy.append(v.copy())
# j = [j1, j2, j3, j4, j5, j6]

dq = []         # qf - q0
vel_max = []    # sqrt(acc_max * dq)
times = []      # contains [tb, T, tf] for each joint

q = []          # Displacement curve
v = []          # Velocity curve
a = []          # Acceleration curve
bq = []         # Buffer Displacement curve
bv = []         # Buffer Velocity curve
ba = []         # Buffer Acceleration curve
time_intervals = []     # Time Intervals from 0 to tf for each Joint to move

Q0 = 0          # Index for q0
QF = 1          # Index for qf
VEL = 2         # Index for velocity max
ACC = 3         # Index for acceleration max
T_B = 4         # Index for tb
T_T = 5         # Index for T
T_F = 6         # Index for tf

flag = 0

dwell_old = 0
new_dwell = 0
dq_real = []

wrist = ['z', 'x', 'z']     # Euler's Wrist Configuration (xyz, xzy, xyx, xzx, yxz, yzx, yxy, yzy, zxy, zyx, zyz, zxz)

def get_times_trajectory(j):
    global dq
    dq = []
    vel_max = []
    times = []
    for index in range(len(j)):
        dq.append(j[index][QF]-j[index][Q0])
        vel_max.append(np.sqrt(j[index][ACC] * dq[index]))

    for index in range(len(vel_max)):
        if vel_max[index] <= j[index][VEL]:    # Triangular Profile
            tb = np.sqrt(dq[index] / j[index][ACC])
            T = tb
            tf = 2 * tb
        elif vel_max[index] > j[index][VEL]:      # Trapezoidal Profile
            tb = j[index][VEL] / j[index][ACC]      # tb = vel_max / acc_max
            T = dq[index] / j[index][VEL]         # T = tf-tb = (qf-q0) / vel
            tf = T + tb                         # tf = T + tb
        else:
            tb, T, tf = None, None, None
        j[index].append(tb)
        j[index].append(T)
        j[index].append(tf)
        times.append([tb, T, tf])

    for index in range(len(j)):
        print(f'For Joint {index+1} : tb = {j[index][4]:.2f},  tf-tb = {j[index][5]:.2f},  tf = {j[index][6]:.2f}')

    return j, times

def plan_trajectories(j, numerical=False):
    t0 = 0
    time_intervals.clear()
    time_stamp = 0
    q.clear()
    v.clear()
    a.clear()
    for index in range(len(j)):
        if numerical is False:
            time_stamp = np.linspace(0, j[index][T_F], 20)
            time_intervals.append(list(time_stamp))
        elif numerical is True:
            time_stamp = np.arange(0, j[index][T_F] + dt, dt)
            time_intervals.append(list(time_stamp))
        for t in time_stamp:
            if j[index][T_B] == j[index][T_T]:  # If it is Triangular Profile
                if t <= j[index][T_B]:
                    qi = j[index][Q0] + (0.5 * j[index][ACC] * (t-t0)**2)
                    q02 = qi
                    vi = j[index][ACC] * t
                    v02 = vi
                    ai = j[index][ACC]

                elif j[index][T_B] < t <= j[index][T_F]:
                    vi = j[index][ACC] * (j[index][T_F] - t)
                    qi = j[index][QF] - (0.5 * j[index][ACC] * (t - j[index][T_F]) ** 2)
                    ai = -j[index][ACC]
            else:                              # If it is Trapezoidal Profile
                if t <= j[index][T_B]:
                    qi = j[index][Q0] + (0.5 * j[index][ACC] * (t-t0)**2)
                    q02 = qi
                    vi = j[index][ACC] * t
                    v02 = vi
                    ai = j[index][ACC]

                elif j[index][T_B] < t <= j[index][T_T]:
                    vi = j[index][VEL]
                    qi = q02 + v02 * (t - j[index][T_B])
                    ai = 0

                elif t > j[index][T_T]:
                    vi = j[index][ACC] * (j[index][T_F] - t)
                    qi = j[index][QF] - (0.5 * j[index][ACC] * (t - j[index][T_F]) ** 2)
                    ai = -j[index][ACC]

            bq.append(qi)
            bv.append(vi)
            ba.append(ai)
        q.append(bq.copy())
        v.append(bv.copy())
        a.append(ba.copy())
        bq.clear()
        bv.clear()
        ba.clear()
    return time_intervals, q, v, a, j

def plot_trajectories(t, q, v, a, j, titlel='Plot'):

    pack = [q, v, a]
    tf = []
    for index in range(len(j)):
        tf.append(j[index][T_F])

    count = 0
    figure(figsize=(26, 8))
    suptitle(titlel)
    for index in range(len(j)):
        for i in range(len(pack)):
            count += 1
            if i == 0:
                label = r"q(t) ($\degree$)"
                ylimit = [0, max(pack[i][index]) + 10]
                slabel = "q"
                mlabel = r"$q_{2}$"
            elif i == 1:
                label = r"v(t) ($\degree$/s)"
                ylimit = [0, max(pack[i][index]) + 1]
                slabel = "v"
                mlabel = r"$v_{max}$"

            elif i == 2:
                label = r"a(t) ($\degree/s^2$)"
                ylimit = [min(pack[i][index])-1, max(pack[i][index])+1]
                slabel = "a"
            subplot(int(3*100+3*10+count))
            plot(t[index], pack[i][index], linewidth=2, label=slabel)  # pack[i][index]
            xlabel('t (s)', fontsize=18)
            ylabel(label, fontsize=18)
            grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
            grid(True)
            xlim([0, max(tf)])
            ylim(ylimit)
            if i == 0:
                hlines(j[index][QF], 0, j[index][T_F], linestyles='--', color='r', label=mlabel)
                vlines([j[index][T_B], j[index][T_T]], 0, ylimit[1], linestyles='--', linewidth=2)
            elif i == 1:
                hlines(j[index][VEL], 0, j[index][T_F], linestyles='--', color='r', label=mlabel)
                vlines([j[index][T_B], j[index][T_T]], 0, ylimit[1], linestyles='--', linewidth=2)
            elif i == 2:
                vlines([j[index][T_B], j[index][T_T]], 0, ylimit, linestyles='--', linewidth=2)
            legend()
            if count == 9:
                show()
                count = 0
                global flag
                flag = flag + 1
                figure(figsize=(26, 8))
                if flag <= 1:
                    suptitle(titlel)

    show()
    flag = 0
    return None

def synchronize_trajectories(j):
    tb_all = []  # tb (Rise) time of all the joints
    T_all = []  # tf-tb time of all the joints
    tf_all = []  # tf (final) time of all the joints
    dwell_all = []  # tf-tb-tb (dwell) time of all the joints
    for index in range(len(j)):
        tb_all.append(j[index][T_B])
        T_all.append(j[index][T_T])
        tf_all.append(j[index][T_F])
        dwell_all.append(T_all[index] - tb_all[index])

    tb = max(tb_all)
    print("")
    for index in range(len(tb_all)):
        print(f"Joint {index+1} rise time : {tb_all[index]:.2f}")
    print(f"Synced rise time : {tb:.2f}\n")

    dwell = max(dwell_all)
    for index in range(len(dwell_all)):
        print(f"Joint {index+1} dwell time : {dwell_all[index]:.2f}")
    print(f"Synced dwell time : {dwell:.2f}\n")

    T = dwell + tb
    tf = T + tb
    print(f"Synchronized trajectory time = {tb:.2f} + {dwell:.2f} + {tb:.2f} = {tf:.2f} sec\n")

    # Update the new tb, T, tf, vel_max, acc_max values to all the Joints in j variable
    for index in range(len(j)):
        j[index][T_B] = tb
        j[index][T_T] = T
        j[index][T_F] = tf
        j[index][VEL] = (j[index][QF] - j[index][Q0]) / T
        j[index][ACC] = j[index][VEL] / tb

    for i in range(len(j)):
        print(f"Joint {i+1} velocity modified from {j_copy[i][VEL]:.2f} to {j[i][VEL]:.2f} and acceleration from {j_copy[i][ACC]:.2f} to {j[i][ACC]:.2f}")

    return j

def calc_propagation_error(j):
    global dq_real
    dq_real = []
    for index in range(len(j)):
        if (j[index][T_B] * 1000) % (dt * 1000) != 0:
            global dwell_old, new_dwell
            dwell_old = (j[index][T_T] - j[index][T_B])
            # print(j[1][T_B] * 1000, dt * 1000)
            j[index][T_B] = (((dt * 1000) - ((j[index][T_B] * 1000) % (dt * 1000))) + (j[index][T_B] * 1000)) / 1000
            new_dwell = (((dt * 1000) - ((dwell_old * 1000) % (dt * 1000))) + (dwell_old * 1000)) / 1000
            j[index][T_T] = new_dwell + j[index][T_B]
            j[index][T_F] = j[index][T_T] + j[index][T_B]
            j[index][VEL] = j[index][ACC] * j[index][T_B]
            # jn[index][VEL] = (jn[index][QF] - jn[index][Q0]) / (jn[index][T_T])
            # jn[index][ACC] = jn[index][VEL] / jn[index][T_B]
            dq_real.append(float("{:.2f}".format(j[index][VEL] * (j[index][T_F] - j[index][T_B]))))

        else:
            # if the tb, T, tf are divisible by dt then
            # we don't need to update new values of tb, T, tf
            # Because tb, T, tf will remain same
            return j


    print(f"The new values for each joint at Controller Freq = {1/dt}Hz")
    for index in range(len(j)):
        print(f"For Joint {index+1} new tb = {j[index][T_B]:.2f}, new tf-tb = {j[index][T_T]:.2f}, new tf = {j[index][T_F]:.2f}, new dq = {dq_real[index]:.2f}, new vel_max = {j[index][VEL]:.2f}, old acc_max = {j[index][ACC]:.2f}")

    pos_real = forward(dq_real, d2=1, d6=1)
    pos_actual = forward(dq, d2=1, d6=1)

    print("\ndq_real is : ", dq_real)
    print("dq_actual is : ", dq)
    print("Position real is : ", pos_real)
    print("Position actual is : ", pos_actual)

    propagated_error = []
    for pos in range(len(pos_real)):
        propagated_error.append(float("{:.2f}".format(pos_real[pos] - pos_actual[pos])))

    print(f"\nPropagated Error in End effector for Controller Freq {1/dt}Hz is : {propagated_error}")

    return j

def forward(joint_params, d2, d6):
    """Forward Kinematics Solution For Stanford Manipulator
    with Spherical Wrist in any valid euler_wrist arrangement
    :param joint_params: Contains t1, t2, d3, t4, t5, t6 and d3(d3 is not an angle its displacement)
    :param d2: Link 2 length
    :param d6: Link 6 length"""

    t1, t2, d3, t4, t5, t6 = joint_params
    t1, t2, d3, t4, t5, t6 = np.deg2rad(t1), np.deg2rad(t2), d3, np.deg2rad(t4), np.deg2rad(t5), np.deg2rad(t6)

    t01 = rotz(t1)                                  # Rotation around Z axis
    t12 = rotx(t2) * tranx(d2)                      # Rotation around Z axis and Translation along X axis
    t23 = tranz(d3)                                 # Translation along Z axis
    t34 = rotz(t4)                                  # Rotation around Z axis    #rotz(t4)
    t45 = rotx(t5)                                  # Rotation around X axis    #rotx(t5)
    t56 = rotz(t6) * tranz(d6)                      # Rotation around Z axis and Translation along Z axis   #rotz(t6)
    t03 = t01 * t12 * t23                           # Transformation from 0 to 3
    t36 = t34 * t45 * t56                           # Transformation from 3 to 6
    t06 = t03 * t36                                 # Transformation from 0 to 6 (Complete Forward Kinematic Solution)
    X = t06[0, 3]                                   # X position of end effector
    Y = t06[1, 3]                                   # Y position of end effector
    Z = t06[2, 3]                                   # Z position of end effector
    return [float("{:.2f}".format(X)), float("{:.2f}".format(Y)), float("{:.2f}".format(Z))]

def plan_plot_polynomial_trajectories(j_1, tif, titlel='Plot'):
    ax = [[0.0]]*len(j_1)
    for joint in range(len(j_1)):
        Q = Matrix([[j_1[joint][0]],
                    [j_1[joint][1]],
                    [j_1[joint][2]],
                    [j_1[joint][3]],
                    [j_1[joint][4]],
                    [j_1[joint][5]]])
        t0 = tif[joint][0]
        tf = tif[joint][1]
        M = Matrix([[t0 ** 5,       t0 ** 4,      t0 ** 3,     t0 ** 2,     t0 ** 1,    1],
                    [tf ** 5,       tf ** 4,      tf ** 3,     tf ** 2,     tf ** 1,    1],
                    [5 * t0 ** 4,   4 * t0 ** 3,  3 * t0 ** 2, 2 * t0 ** 1, 1,          0],
                    [5 * tf ** 4,   4 * tf ** 3,  3 * tf ** 2, 2 * tf ** 1, 1,          0],
                    [20 * t0 ** 3,  12 * t0 ** 2, 6 * t0 ** 1, 2,           0,          0],
                    [20 * tf ** 3,  12 * tf ** 2, 6 * tf ** 1, 2,           0,          0]])
        A_vals =  M.inv() * Q
        ax[joint] = A_vals.tolist()
        print(f"For Joint {joint+1} a5, a4, a3, a2, a1, a0 are : {A_vals.tolist()}")

    q.clear()
    v.clear()
    a.clear()
    time_stamp = 0
    time_intervals.clear()
    #print(ax[0][0][0])
    for joint in range(len(j_1)):
        time_stamp = np.linspace(0, tif[joint][1], 20)
        time_intervals.append(list(time_stamp))
        for t in time_stamp:
            qi = ax[joint][0][0]*t**5 + ax[joint][1][0]*t**4 + ax[joint][2][0]*t**3 + ax[joint][3][0]*t**2 + ax[joint][4][0]*t + ax[joint][5][0]
            vi = 5*ax[joint][0][0]*t**4 + 4*ax[joint][1][0]*t**3 + 3*ax[joint][2][0]*t**2 + 2*ax[joint][3][0]*t + ax[joint][4][0]
            ai = 20*ax[joint][0][0]*t**3 + 12*ax[joint][1][0]*t**2 + 6*ax[joint][2][0]*t + 2*ax[joint][3][0]
            bq.append(qi)
            bv.append(vi)
            ba.append(ai)
        q.append(bq.copy())
        v.append(bv.copy())
        a.append(ba.copy())
        bq.clear()
        bv.clear()
        ba.clear()

    for joint in range(len(j_1)):
        figure(figsize=(26,8))
        suptitle(titlel)
        subplot(132)
        plot(time_intervals[joint], v[joint], linewidth=2, label="v")
        xlabel('t (s)', fontsize=18)
        ylabel(r'v(t) ($\degree$/s)', fontsize=18)
        grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        grid(True)
        legend()

        subplot(131)
        plot(time_intervals[joint], q[joint], linewidth=2, label="q")
        xlabel('t (s)', fontsize=18)
        ylabel(r'q(t) ($\degree$)', fontsize=18)
        grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        grid(True)
        legend()

        subplot(133)
        plot(time_intervals[joint], a[joint], linewidth=3, label="a")
        xlabel('t (s)', fontsize=18)
        ylabel(r'a(t) ($\degree^2$/s)', fontsize=18)
        grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        grid(True)
        legend()
        print(f"\n\nFor Joint {joint+1} Plots of Position, Velocity and Acceleration : ")
        show()
    return None

def trajectories_4_consecutive_points(j_1, tif, titlel='Plot'):
    ax = [[0.0]] * len(j_1)
    for joint in range(len(j_1)):
        Q = Matrix([[j_1[joint][0]],
                    [j_1[joint][1]],
                    [j_1[joint][1]],
                    [j_1[joint][2]],
                    [j_1[joint][2]],
                    [j_1[joint][3]],
                    [j_1[joint][4]],
                    [j_1[joint][5]],
                    [0],
                    [0],
                    [0],
                    [0]])
        t0 = tif[joint][0]
        t1 = tif[joint][1]
        t2 = tif[joint][2]
        tf = tif[joint][3]

        M = Matrix([[t0 ** 3, t0 ** 2, t0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [t1 ** 3, t1 ** 2, t1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, t1 ** 3, t1 ** 2, t1, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, t2 ** 3, t2 ** 2, t2, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, t2 ** 3, t2 ** 2, t2, 1],
                    [0, 0, 0, 0, 0, 0, 0, 0, tf ** 3, tf ** 2, tf, 1],
                    [3 * t0 ** 2, 2 * t0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 3 * tf ** 2, 2 * tf, 1, 0],
                    [3 * t1 ** 2, 2 * t1, 1, 0, -3 * t1 ** 2, -2 * t1, -1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 3 * t2 ** 2, 2 * t2, 1, 0, -3 * t2 ** 2, -2 * t2, -1, 0],
                    [6 * t1, 2, 0, 0, -6 * t1, -2, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 6 * t2, 2, 0, 0, -6 * t2, -2, 0, 0]])

        A_vals = M.inv() * Q
        ax[joint] = A_vals.tolist()
        print(f"For Joint {joint + 1} a31, a21, a11, a01, a32, a22, a12, a02, a33, a23, a13, a03 are : "
              f"{A_vals.tolist()}")

    q.clear()
    v.clear()
    a.clear()
    time_stamp = 0
    time_intervals.clear()
    #print(ax[0][0][0])
    for joint in range(len(j_1)):
        time_stamp = np.linspace(0, tif[joint][3], 3000)
        time_intervals.append(list(time_stamp))
        for t in time_stamp:
            if t <= tif[joint][1]:
                qi = ax[joint][0][0]*t**3 + ax[joint][1][0]*t**2 + ax[joint][2][0]*t + ax[joint][3][0]
                vi = 3*ax[joint][0][0]*t**2 + ax[joint][1][0]*t + ax[joint][2][0]
                ai = 6*ax[joint][0][0]*t + ax[joint][1][0]
            elif tif[joint][1] < t <= tif[joint][2]:
                qi = ax[joint][4][0]*t**3 + ax[joint][5][0]*t**2 + ax[joint][6][0]*t + ax[joint][7][0]
                vi = 3*ax[joint][4][0]*t**2 + ax[joint][5][0]*t + ax[joint][6][0]
                ai = 6*ax[joint][4][0]*t + ax[joint][5][0]
            elif tif[joint][2] < t <= tif[joint][3]:
                qi = ax[joint][8][0]*t**3 + ax[joint][9][0]*t**2 + ax[joint][10][0]*t + ax[joint][11][0]
                vi = 3*ax[joint][8][0]*t**2 + ax[joint][9][0]*t + ax[joint][10][0]
                ai = 6*ax[joint][8][0]*t + ax[joint][9][0]

            bq.append(qi)
            bv.append(vi)
            ba.append(ai)
        q.append(bq.copy())
        v.append(bv.copy())
        a.append(ba.copy())
        bq.clear()
        bv.clear()
        ba.clear()

    for joint in range(len(j_1)):
        figure(figsize=(26,8))
        suptitle(titlel)
        subplot(132)
        plot(time_intervals[joint], v[joint], linewidth=2, label="v")
        xlabel('t (s)', fontsize=18)
        ylabel(r'v(t) ($\degree$/s)', fontsize=18)
        grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        grid(True)
        legend()

        subplot(131)
        plot(time_intervals[joint], q[joint], linewidth=2, label="q")
        xlabel('t (s)', fontsize=18)
        ylabel(r'q(t) ($\degree$)', fontsize=18)
        grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        grid(True)
        legend()

        subplot(133)
        plot(time_intervals[joint], a[joint], linewidth=3, label="a")
        xlabel('t (s)', fontsize=18)
        ylabel(r'a(t) ($\degree^2$/s)', fontsize=18)
        grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        grid(True)
        legend()
        print(f"\n\nFor Joint {joint+1} q0:{j_1[joint][0]}, q1:{j_1[joint][1]}, q1:{j_1[joint][2]}, q1:{j_1[joint][3]}")
        show()

    return None

def rotz(t):
    """Rotates the frame about Z Axis"""
    rot = Matrix([[cos(t), -sin(t), 0, 0],
                  [sin(t), cos(t),  0, 0],
                  [0,       0,      1, 0],
                  [0,       0,      0, 1]])
    return rot

def rotx(t):
    """Rotates the frame about X Axis"""
    rot = Matrix([[1,       0,      0, 0],
                  [0, cos(t), -sin(t), 0],
                  [0, sin(t), cos(t),  0],
                  [0,   0,      0,     1]])
    return rot

def roty(t):
    """Rotates the frame about Y Axis"""
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

if __name__ == '__main__':

    start_time = time.time()
    tprint("Task   -   1")
    j, times = get_times_trajectory(j)
    time_intervals, q, v, a, j = plan_trajectories(j)
    plot_trajectories(time_intervals, q, v, a, j, titlel='Task-1 Trajectories without Synchronization for Stanford '
                                                         'Robot model')
    print("-----------------END of Task 1-----------------\n\n")

    tprint("Task   -   2")
    j = synchronize_trajectories(j)
    time_intervals, q, v, a, j = plan_trajectories(j)
    plot_trajectories(time_intervals, q, v, a, j, titlel='Task-2 Synchronized Trajectories')
    print("-----------------END of Task 2-----------------\n\n")

    tprint("Task   -   3")
    j = calc_propagation_error(j)
    time_intervals, q, v, a, j = plan_trajectories(j, numerical=True)
    plot_trajectories(time_intervals, q, v, a, j, titlel='Task-3 Propagated Error Calculation in the End Effector')
    print("-----------------END of Task 3-----------------\n\n")

    tprint("Task   -   4")
    j, times = get_times_trajectory(j)
    j = synchronize_trajectories(j)
    j = calc_propagation_error(j)
    time_intervals, q, v, a, j = plan_trajectories(j, numerical=True)
    plot_trajectories(time_intervals, q, v, a, j, titlel='Task-4 Synchronized Trajectories for Numerical Control')
    print("-----------------END of Task 4-----------------\n\n")

    tprint("Task   -   5")
    j1 = [0, 20, 0, 0.1, 0, 2]      # Joint 1 = [q0, qf, v0, vf, a0, af]
    j2 = [0, 30, 0, 0.1, 0, 2]      # Joint 2 = [q0, qf, v0, vf, a0, af]
    j3 = [0, 40, 0, 0.1, 0, 2]      # Joint 3 = [q0, qf, v0, vf, a0, af]
    j4 = [0, 50, 0, 0.1, 0, 2]      # Joint 4 = [q0, qf, v0, vf, a0, af]
    j5 = [0, 60, 0, 0.1, 0, 2]      # Joint 5 = [q0, qf, v0, vf, a0, af]
    j6 = [0, 70, 0, 0.1, 0, 2]      # Joint 6 = [q0, qf, v0, vf, a0, af]
    j_1 = [j1, j2, j3, j4, j5, j6]
    tif = [[0, 5]] * 6
    plan_plot_polynomial_trajectories(j_1, tif, titlel='Task-5 Polynomial Trajectory for 2 points for 6 Joints')
    print("-----------------END of Task 5-----------------\n\n")

    tprint("Task   -   6")
    j1 = [0, 10, 20, 40, 0, 5]      # Joint 1 = [q0, q1, q2, qf, v0, vf]
    j2 = [0, 10, 20, 40, 0, 0]      # Joint 2 = [q0, q1, q2, qf, v0, vf]
    j3 = [0, 10, 20, 40, 0, 0]      # Joint 3 = [q0, q1, q2, qf, v0, vf]
    j4 = [0, 10, 20, 40, 0, 0]      # Joint 4 = [q0, q1, q2, qf, v0, vf]
    j5 = [0, 10, 20, 40, 0, 0]      # Joint 5 = [q0, q1, q2, qf, v0, vf]
    j6 = [0, 10, 20, 40, 0, 0]      # Joint 6 = [q0, q1, q2, qf, v0, vf]
    j_1 = [j1, j2, j3, j4, j5, j6]
    t0 = 0; t1 = 5; t2 = 10; tf = 15
    tif = [[t0, t1, t2, tf]] * 6
    trajectories_4_consecutive_points(j_1, tif, titlel='Task-6 Polynomial Trajectory for 4 points for 6 Joints')
    print("-----------------END of Task 6-----------------\n\n")

    print(f"Executed this program in : {time.time() - start_time} seconds\n")
