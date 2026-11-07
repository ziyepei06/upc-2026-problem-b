"""Quasi-steady blow-down of a rigid, well-mixed module through a sharp orifice.

The model is Eqs. (1)-(3) of the paper: mass balance for the reservoir,
choked and subsonic orifice speed, and a first-order wall heat transfer.
Integrated with explicit RK4 and a fixed step small enough that halving it
changes t_0.3 by less than 0.1 % (see Section 5).

Run:  python3 blowdown.py            # nominal case, prints t_0.3
"""

import math

P_ATM = 101325.0
GAMMA = 1.4
R_AIR = 287.05


def exit_speed(p, p_a, T):
    """Exit speed (m/s) from the reservoir state, choked or subsonic."""
    ratio = p_a / p if p > 0 else 1.0
    critical = (2.0 / (GAMMA + 1.0)) ** (GAMMA / (GAMMA - 1.0))
    if ratio <= critical:
        return math.sqrt(GAMMA * R_AIR * T) * math.sqrt(
            2.0 / (GAMMA + 1.0))
    return math.sqrt(2.0 * GAMMA / (GAMMA - 1.0) * R_AIR * T
                     * (1.0 - ratio ** ((GAMMA - 1.0) / GAMMA)))


def rhs(p, T, A, V, Cd, p_a, T_w, tau_w, dt):
    """dp/dt and dT/dt for the reservoir."""
    u = exit_speed(p, p_a, T)
    dp = -GAMMA * p * A * Cd / V * u
    dT = (GAMMA - 1.0) / GAMMA * T / p * dp - (T - T_w) / tau_w
    return dp, dT


def blowdown(p0=P_ATM, V=628.0, d=0.10, Cd=0.82, T0=293.0, T_w=293.0,
             tau_w=300.0, dt=0.5):
    """Return (t, p, T) arrays until p falls to 0.3 p0."""
    A = math.pi * d ** 2 / 4.0
    p, T, t = p0, T0, 0.0
    times, pressures, temps = [t], [p], [T]
    while p > 0.3 * p0:
        k1p, k1T = rhs(p, T, A, V, Cd, 0.0, T_w, tau_w, dt)
        k2p, k2T = rhs(p + dt / 2 * k1p, T + dt / 2 * k1T, A, V, Cd, 0.0,
                       T_w, tau_w, dt)
        k3p, k3T = rhs(p + dt / 2 * k2p, T + dt / 2 * k2T, A, V, Cd, 0.0,
                       T_w, tau_w, dt)
        k4p, k4T = rhs(p + dt * k3p, T + dt * k3T, A, V, Cd, 0.0, T_w,
                       tau_w, dt)
        p += dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        T += dt / 6 * (k1T + 2 * k2T + 2 * k3T + k4T)
        t += dt
        times.append(t)
        pressures.append(p)
        temps.append(T)
    return times, pressures, temps


if __name__ == "__main__":
    t, p, T = blowdown()
    print(f"t_0.3 = {t[-1]:.1f} s for d = 10 cm, V = 628 m^3")
