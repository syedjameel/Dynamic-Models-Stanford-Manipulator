from trajectory_planning import *
from SimpleTranformations import setup_symbolic

def plan_plot_trajectories(j_1, tif, titlel='Plot'):
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

def forward2(joint_params, d2, d6, symbolic):
    """Forward Kinematics Solution For Stanford Manipulator
    with Spherical Wrist in any valid euler_wrist arrangement
    :param joint_params: Contains t1, t2, d3, t4, t5, t6 and d3(d3 is not an angle its displacement)
    :param d2: Link 2 length
    :param d6: Link 6 length"""
    if symbolic is False:
        t1, t2, t3, t4, t5, t6 = joint_params
        d2 = d2; d6 = d6
        from SimpleTranformations import rotz, rotx, tranx, tranz

    if symbolic is True:
        from sympy.solvers import solve
        setup_symbolic(symbolic)
        from SimpleTranformations import rotz, rotx, roty, tranx, tranz
        from sympy import symbols
        t1 = symbols('t1')
        t2 = symbols('t2')
        t4 = symbols('t4')
        t5 = symbols('t5')
        t6 = symbols('t6')
        d2 = symbols('d2')
        d3 = symbols('d3')
        d6 = symbols('d6')

    t01 = rotz(t1)               # Rotation around Z axis
    t12 = rotx(t2) * tranx(d2)   # Rotation around Z axis and Translation along X axis
    t23 = tranz(d3)              # Translation along Z axis
    t34 = rotz(t4)               # Rotation around Z axis    #rotz(t4)
    t45 = rotx(t5)               # Rotation around X axis    #rotx(t5)
    t56 = rotz(t6) * tranz(d6)   # Rotation around Z axis and Translation along Z axis
    t36 = t34 * t45 * t56         # Transformation from 3 to 6
    t02 = t01 * t12
    t03 = t02 * t23
    t04 = t03 * t34
    t05 = t04 * t45
    t06 = t05 * t56              # Transformation from 0 to 6
    return t01, t12, t23, t34, t45, t56

