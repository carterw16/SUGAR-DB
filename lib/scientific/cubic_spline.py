import numpy as np
from scipy.interpolate import CubicSpline


def create_cubicSpline(Xbounds, Ybounds):
    # create cubic spline given a set of bounds and return the cubic spline coefficients and the spline
    cs = CubicSpline(Xbounds, Ybounds, bc_type='natural')
    a = np.array([cs.c[0, 0], cs.c[0, 1], cs.c[0, 2]])
    b = np.array([cs.c[1, 0], cs.c[1, 1], cs.c[1, 2]])
    c = np.array([cs.c[2, 0], cs.c[2, 1], cs.c[2, 2]])
    d = np.array([cs.c[3, 0], cs.c[3, 1], cs.c[3, 2]])
    return cs, a, b, c, d


def cubicSpline_partials(a, b, c, Xr, Xi, SF, p):
    X = np.abs(Xr + Xi * 1j)
    dX_dXr = Xr / X
    dX_dXi = Xi / X
    dF_dX = -(X - p) * (3 * a * (X - p) + (2 * b)) - c
    # return the partial derivatives defined by the spline region
    dF_dXr = dF_dX * dX_dXr
    dF_dXi = dF_dX * dX_dXi
    return dF_dXr, dF_dXi

# def cubicSpline_partials(a, b, c, X, p):
#     dF_dX = (-3 * a * (X - p)**2) - (2 * b * (X - p)) - c
#     # return the partial derivatives defined by the spline region
#     # dF_dXr = dF_dX * dX_dXr
#     # dF_dXi = dF_dX * dX_dXi
#     return dF_dX


def linear_partials(m, Xr, Xi, SF):
    dF_dXr = -SF * m * Xr
    dF_dXi = -SF * m * Xi
    return dF_dXr, dF_dXi
