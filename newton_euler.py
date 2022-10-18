#########################
# Newton Euler Formulas #
#########################

def w_i(joint_type, w_i_1, R_i_1, dq_i, u_i_1):
    """Calculate the Angular Velocity on link i"""
    if joint_type == 'p':
        w = R_i_1.T * w_i_1
    elif joint_type == 'r':
        w = R_i_1.T * (w_i_1 + dq_i * u_i_1)

    return w

def dw_i(joint_type, w_i_1, dw_i_1, R_i_1, dq_i, ddq_i, u_i_1):
    """Calculates the Angular Acceleration on Link i"""
    if joint_type == 'p':
        dw = R_i_1.T * dw_i_1
    elif joint_type == 'r':
        dw = R_i_1.T * (dw_i_1 + (ddq_i * u_i_1) + (dq_i * w_i_1).cross(u_i_1))

    return dw

def ddp_i(joint_type, R_i_1, ddq_i, ddp_i_1, dq_i, w_i, dw_i, r_i_1_1, u_i_1):
    """Calculates the Linear Acceleration"""
    if joint_type == 'p':
        ddp = R_i_1.T *(ddp_i_1 + ddq_i*u_i_1) + (2*dq_i*w_i).cross((R_i_1.T*u_i_1)) \
                + dw_i.cross(r_i_1_1) + w_i.cross(w_i.cross(r_i_1_1))
    elif joint_type == 'r':
        ddp = R_i_1.T *ddp_i_1 + dw_i.cross(r_i_1_1) + w_i.cross(w_i.cross(r_i_1_1))

    return ddp

def ddp_ci(ddp_i, dw_i, r_i_ci, w_i):
    """Calculates the Linear Acceleration at Centre of Mass"""
    ddp_c = ddp_i + dw_i.cross(r_i_ci) + w_i.cross(w_i.cross(r_i_ci))

    return ddp_c

def forces_moments(R_i1, f_i1, m_i, ddp_ci, meu_i1, r_i_1_1, r_i_ci, dw_i, I_i, w_i):
    """Calculates the Forces and Torques on Link i"""
    f_i = R_i1*f_i1 + m_i*ddp_ci
    meu_i = R_i1*meu_i1 - f_i.cross((r_i_1_1+r_i_ci)) \
            + (R_i1*f_i1).cross(r_i_ci) + I_i*dw_i + w_i.cross((I_i*w_i))

    return f_i, meu_i

def torques_joint(joint_type, f_i, R_i, u_i_1, meu_i):
    if joint_type == 'p':
        tau_i = f_i.T * R_i.T * u_i_1
    elif joint_type == 'r':
        tau_i = meu_i.T * R_i.T * u_i_1

    return tau_i


