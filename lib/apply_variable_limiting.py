import numpy as np


def find_max_min(Vsol, Vnom, Nodes):
    V = []

    V = abs(Vsol[Nodes.Vr_index] + 1j * Vsol[Nodes.Vi_index])
    Vr_mag = np.array(Nodes.Vr_mag, dtype=float, ndmin=2).T

    Vsol_pu = np.divide(V, Vr_mag)

    # find the maximum node voltage and its location
    Vmax = np.amax(Vsol_pu)
    loc_max = np.argmax(Vsol_pu)

    # find the minimum node voltage and its location
    Vmin = np.amin(Vsol_pu)
    loc_min = np.argmin(Vsol_pu)

    return Vmax, Vmin, loc_max, loc_min


def apply_variable_limiting(Vsol, Vmax_lim, Vmin_lim, sigma, sigma_count,
                            errStore, nodes, Nodes):
    sigma_max = 1

    # update nominal voltages/ setpoints
    Vnom = []

    for node in nodes:
        if not node.isTriplex:
            node.update_Vnom(Vsol)
        Vnom.append(node.Vnom)

    # find maximum and minimum voltages in the solution vector
    Vmax, Vmin, loc_max, loc_min = find_max_min(Vsol, Vnom, Nodes)

    # check if the step is too large
    if Vmax > Vmax_lim and sigma > 1e-2:
        sigma *= 0.5
    elif Vmin < Vmin_lim and sigma > 1e-2:
        sigma *= 0.5

    sigma_count += 1

    # if variable limiting applied at least once, check if error is monotonically decreasing
    # if monotonically decreasing and sigma is greater than 1 then scale back up
    if errStore[-1] < errStore[-2]:
        sigma_count += 1
    else:
        sigma_count = 0

    if sigma_count > 5:
        sigma /= 0.75
        min(sigma, sigma_max)
