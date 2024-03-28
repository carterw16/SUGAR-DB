from types import SimpleNamespace

import numpy as np
import scipy.linalg as sl

# A Python port of the MATLAB Polyfix Function by Are Mjaavatten, Telemark University College, Norway, November 2015


def polyfix(x,
            y,
            n,
            xfix=np.array([]),
            yfix=np.array([]),
            xder=None,
            dydx=None):
    # ensure that x and y inputs are column vectors
    x = np.reshape(x, (-1, 1))
    y = np.reshape(y, (-1, 1))

    nfit = x.size

    if not y.size == nfit:
        raise ValueError("x and y must be the same size")

    nfix = xfix.size

    if not yfix.size == nfix:
        raise ValueError("xfix and yfix must be the same size")

    if xder is None:
        xder = np.array([])
        dydx = np.array([])
    else:
        xder = np.reshape(xder, (-1, 1))
        dydx = np.reshape(dydx, (-1, 1))

    nder = xder.size

    if not dydx.size == nder:
        raise ValueError("xder and dydx must be the same size")

    nspec = nfix + nder

    try:
        specval = np.concatenate((yfix, dydx.flatten()), axis=1)
    except np.AxisError:
        specval = np.concatenate((yfix, dydx.flatten()), axis=0).reshape(
            (-1, 1))

    # First find A and pc such that A*pc = specval
    A = np.zeros((nspec, n + 1))

    # Specific values of y
    if not (nfix == 0):
        for i in range(n + 1):
            A[0:nfix, i] = (np.ones((nfix, 1)) * xfix**((n - i))).T

    # Specified values of dydx
    if nder > 0:
        for i in range(n):
            A[nfix:nfix + nder, i] = ((n - i + 1) * np.ones(
                (nder, 1)) * xder**((n - i))).T

    if nfix > 0:
        lastcol = n + 1
        nmin = nspec - 1
    else:
        lastcol = n  # if only derivatives p(n+1) is arbitrary
        nmin = nspec

    if n < nmin:
        raise ValueError(
            "Polynmoial degree too low. Cannot match all constraints")

    # Find the unique polynomial of degree nmin that fits the constraints

    firstcol = n - nmin
    pc0 = sl.solve(A[:, firstcol:lastcol], specval)

    # Extend and pad

    pc = np.zeros((n + 1, 1))

    pc[firstcol:lastcol] = pc0.reshape((-1, 1))

    # Create the Vandermode Matrix
    x_v = x.flatten()
    V = np.vander(x_v, n + 1)

    yfit = y - np.polyval(pc, x)

    # Find the p0 that minimizes (V*p0 - yfit)'*(V*p0-yfit)

    B = sl.null_space(A)
    z = np.linalg.lstsq(np.matmul(V, B), yfit, rcond=None)[0]
    p0 = np.matmul(B, z)
    p = p0.T + pc.T

    # Error Estimation

    R = sl.qr(V, mode='r')
    r = y - np.matmul(V, p.T)
    df = np.max(np.array([0, y.size - (n + 1)]))
    S = {'R': R, 'df': df, 'normr': sl.norm(r)}
    S = SimpleNamespace(**S)

    return p
