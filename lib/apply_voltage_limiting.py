import numpy as np


def apply_voltage_limiting(Vsol, V, Nodes, maxStep__, minStep__, Vr_norm,
                           Vi_norm, VrLimit__, ViLimit__, stamp_dual = False):
    # Apply voltage limiting module here
    # deltaTap = Vsol[Nodes.Tap_index] - V[Nodes.Tap_index]
    deltaVr = Vsol[Nodes.Vr_index] - V[Nodes.Vr_index]
    deltaVi = Vsol[Nodes.Vi_index] - V[Nodes.Vi_index]
    # Convert the norm vectors into array
    if Vr_norm is None:
        Vr_norm = np.array(Nodes.Vr_mag, dtype=float, ndmin=2).T
    if Vi_norm is None:
        Vi_norm = np.array(Nodes.Vi_mag, dtype=float, ndmin=2).T

    # Normalize the solution deltaVr and deltaVi vectors
    delta_Vr_norm = np.divide(deltaVr, Vr_norm)
    delta_Vi_norm = np.divide(deltaVi, Vi_norm)

    # Min-max of the deltaVr and deltaVi
    delta_Vr_norm = np.minimum(delta_Vr_norm, maxStep__)
    delta_Vr_norm = np.maximum(delta_Vr_norm, minStep__)
    delta_Vi_norm = np.minimum(delta_Vi_norm, maxStep__)
    delta_Vi_norm = np.maximum(delta_Vi_norm, minStep__)

    # Convert back to deltaVr
    delta_Vr = np.multiply(delta_Vr_norm, Vr_norm)
    delta_Vi = np.multiply(delta_Vi_norm, Vi_norm)

    # Update V
    V[Nodes.Vr_index] += delta_Vr
    V[Nodes.Vi_index] += delta_Vi
    V[Nodes.Tap_index] = np.copy(Vsol[Nodes.Tap_index])
    # V[Nodes.I_index] = np.copy(Vsol[Nodes.I_index])
    V[Nodes.Q_index] = np.copy(Vsol[Nodes.Q_index])

    if stamp_dual:
        dual_max_step = 3*maxStep__
        dual_min_step = 3*minStep__
        deltaL = Vsol[Nodes.L_index] - V[Nodes.L_index]
        np.fmin(deltaL, dual_min_step * np.ones((np.shape(deltaL)[0], 1)), out=deltaL)
        np.fmax(deltaL, dual_max_step * np.ones((np.shape(deltaL)[0], 1)), out=deltaL)    
        V[Nodes.L_index] += deltaL

    return V
# Q Limiting
# maxQStep__ = 0.3
# minQStep__ = -maxQStep__
# deltaQ = Vsol[Nodes.Q_index] - V[Nodes.Q_index]
# maxQStep = maxQStep__ * np.ones((deltaQ.size, 1))
# minQStep = minQStep__ * np.ones((deltaQ.size, 1))
# np.fmin(deltaQ, maxQStep, out = deltaQ)
# np.fmax(deltaQ, minQStep, out = deltaQ)
# V[Nodes.Q_index] += deltaQ
