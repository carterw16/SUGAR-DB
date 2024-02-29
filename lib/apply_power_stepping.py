import numpy as np


def apply_power_stepping(NR_count, pstepHardLimit, errMax_Store, pstep_count,
                         pstep_oscillates, err_max, tol, pstep__, min_pstep,
                         V_init, V, power_steppingV, lambdas, not_converged,
                         pstep_Store, pstep_mod):
    flag_ps_continue = False
    flag_ps_break = False
    applyPowerStep = False

    try:
        # Apply power stepping for Hard Limit
        if NR_count >= pstepHardLimit:
            applyPowerStep = True

        else:
            # Apply power stepping if divergence occurs or oscillation occurs
            if errMax_Store[-1] > errMax_Store[-2]:
                pstep_count += 1
                pstep_oscillates += 1
            else:
                pstep_count = 0
                pstep_oscillates -= 1

            if pstep_count > 4:
                applyPowerStep = True
                pstep_count = 0

        # Apply Power Stepping
        if applyPowerStep and err_max > tol:
            if pstep__ < min_pstep:
                print("Power Stepping Failed - Case did not converge")
                raise StopIteration  # to break
            else:
                pstep__ *= 0.5
                V = np.copy(V_init)
                if power_steppingV:
                    V = np.multiply(V, np.sqrt(pstep__))
                lambdas = 1
                NR_count = 0
                print("Power Stepping Step Applied: %s" % str(pstep__))
                raise ValueError  # to continue

        if err_max < tol:
            if pstep__ == 1:
                not_converged = False
            elif pstep__ < 1 and not power_steppingV:
                if pstep__ in pstep_Store:
                    pstep_mod = 1 + 0.5 * (pstep_mod - 1)
                else:
                    pstep_Store.append(pstep__)
                # Update pstep__
                pstep__ = pstep__ * pstep_mod
                pstep__ = min(1, pstep__)
                NR_count = 0
            elif pstep__ < 1 and power_steppingV:
                pstep__ *= 1.1
                pstep__ = min(1, pstep__)
                NR_count = 0
    except ValueError:
        flag_ps_continue = True
        flag_ps_break = True

    return (pstep__, flag_ps_continue, flag_ps_break, pstep_mod, V, lambdas,
            NR_count, not_converged, pstep_count, pstep_oscillates, pstep_Store)
