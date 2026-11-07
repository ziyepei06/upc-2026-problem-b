# Canonical numbers for 230B (copy exactly; do not re-round differently)

All positions in metres, times in seconds, angles in radians unless stated.
Scenario: penalty mark at origin, kick at t=0 from ball centre (0,0,0.11);
goal plane x=11.0, opening |y|<=3.66, 0<=z<=2.44; ball radius 0.110 m.

OBSERVATIONS (sigma_pos = 0.010 m per coordinate assumed unless stated)
  r1 = (1.50, 0.45, 0.44) at t1 = 0.05 s
  r2 = (2.96, 0.86, 0.74) at t2 = 0.10 s

## Model
m dv/dt = -0.5 rho A C_D |v| v + 0.5 rho A C_L R (omega x v)/|v| ... (implemented as km (omega x v))
  rho = 1.20 kg/m3, A = pi R^2, R = 0.110 m, m = 0.430 kg
  kd = 0.5 rho A C_D / m, km = 0.5 rho A C_L R / m   (km uses C_L = 1.0, so |omega| is in rad/s)
  spin decay: |omega(t)| = |omega0| exp(-t/tau), tau = 6 s (axis fixed)
  gravity 9.81 m/s2
At v = 31.45 m/s: |a_drag| = 9.9 m/s2 for C_D = 0.19; Re = 4.4e5.

## E1: model-free constant acceleration (6 parameters: v0, a)
  a = 2[(r1-r0)t2 - (r2-r0)t1] / (t1 t2 (t1-t2)),  v0 = (r1-r0)/t1 - a t1/2,  r0 = (0,0,0.11)
  a = (-16.0, -16.0, -12.0) m/s2 exactly; v0 = (30.4, 9.4, 6.9) m/s (round to these)
  entry: t = 0.4050 s, y = 2.495, z = 1.920
  sigma (sigma_pos = 10 mm, MC 40000): t 0.028, y 0.594, z 0.599; delta-method 0.024/0.576/0.574
  68% y [1.92, 3.09], z [1.35, 2.51]; 95% y [1.55, 3.50], z [0.95, 2.91]
  P(clean goal) = 0.725, P(wide) = 0.029, P(over bar) = 0.190, P(under) = 0.001

## E2: physics-constrained Bayesian (7 parameters: v0, ln C_D, omega)
  priors: ln C_D ~ N(ln 0.19, 0.15^2) [median 0.19, 68% 0.164-0.220];
          omega_i ~ N(0, 35^2) rad/s i.i.d. (marginal |omega| median 33 rad/s ~ 5.3 rev/s);
          C_L = 1.0 fixed, tau = 6 s fixed
  MAP: v0 = (30.080, 8.973, 6.870), ln C_D = -1.633 (C_D = 0.195), omega = (4.82, -2.11, -19.76) rad/s
     entry t = 0.3864 s, y = 3.023, z = 1.945; residuals (mm): -7.2, -9.6, -0.7, +3.7, +4.8, +0.4; chi2 = 1.9 (2 dof)
  MCMC (20000 kept, 4 chains? no: single chain 40000 draws, burn 8000, thin 2, acceptance 0.27, ESS 870-1250):
     posterior mean v0 = (30.073, 8.976, 6.874); C_D = 0.200 +- 0.030;
     |omega| = 52.6 +- 22.1 rad/s = 8.4 +- 3.5 rev/s
  predictive: mean (t 0.387, y 3.029, z 1.945); std (0.004, 0.279, 0.262)
     68% y [2.75, 3.31], z [1.69, 2.20]; 95% y [2.58, 3.50], z [1.52, 2.38]
     P(clean goal) = 0.894, P(wide) = 0.014, P(over bar) = 0.032, P(under) = 0.000

## HEADLINE (equal-weight model average of E1 and E2)
  entry point: x = 11.00 m, y = 2.77 m, z = 1.93 m; t = 0.398 s
  sigma: y 0.54 m, z 0.47 m   (68% y [2.22, 3.25], z [1.54, 2.33];
                                95% y [1.75, 3.50], z [1.15, 2.69])
  P(clean goal) = 0.81; P(wide) = 0.023; P(over the bar) = 0.109; P(clips frame) = 0.058 (samples: 0.809/0.023/0.109/0.058)
  model spread of the y prediction across five defensible models: 2.495 (constant a)
    to 3.217 (drag only / gravity only) -> 0.72 m; z spread 1.919-2.020 (0.10 m)

## Model comparison (sigma_pos = 10 mm)
  M0 gravity only        k=3  chi2 = 6.46  t=0.3706 y=3.217 z=1.957  AIC 12.46
  M1 gravity + drag      k=4  chi2 = 1.52 (C_D=0.379 free) / 2.67 (prior 0.15)  t=0.388 y=3.217 z=1.879-1.919
  M2 full physics        k=7  chi2 = 1.81  t=0.3864 y=3.012 z=1.946  AIC 15.81
  M3 constant accel      k=6  chi2 = 0     t=0.4050 y=2.495 z=1.920  AIC 12.00
  naive ballistic from the measured velocity at t=0.10 s: t=0.3753 y=3.118 z=2.020
  KEY: drag is antiparallel to v, so it cannot change the ratio y/x; y_goal = 11 v0y/v0x = 3.217
  for any C_D. Only a force perpendicular to v (Magnus) or gravity can. Hence the 0.5-0.7 m
  disagreement is entirely a statement about the unobserved lateral (Magnus) force.

## Error budget (share of the variance of the entry-y prediction, sigma_pos = 10 mm)
  E1: y1 74%, y2 25% (x and z coordinates < 1%)
  E2: measurement noise total 35%; omega_z prior 64%; other < 2%
  (so in E1 the lateral acceleration is data-driven; in E2 it is prior-driven)

## Identifiability
  Scaled spectrum of the data Jacobian (rows / sigma, columns = prior-sd units):
    6-parameter model spans 10^4; 7-parameter model spans 10^16
  Fisher eigenvalues of the 7-parameter model (unscaled) span 10^16; adding the prior
  removes only the flat directions (posterior information spans 10^5).
  Compensation: v0y and omega_z are correlated at -0.89; the data fix v0y - c omega_z.
  Moving along the data-equivalent ridge (Delta chi2 = 1) changes the predicted y by ~0.7 m.

## Prior sensitivity (E2, entry y mean +- sd, P(inside), |omega|, C_D)
  base                        3.002 +- 0.266   0.891   54.3   0.195
  ln C_D prior sd 0.30        3.001 +- 0.267   0.892   52.9   0.215
  ln C_D prior sd 0.50        3.012 +- 0.271   0.895   53.4   0.233
  C_D median 0.17             3.004 +- 0.264   0.907   52.8   0.174
  C_D median 0.22             2.997 +- 0.264   0.920   53.8   0.227
  omega prior sd 17 rad/s     3.158 +- 0.153   0.992   26.5   0.194
  omega prior sd 70 rad/s     2.817 +- 0.358   0.797   91.7   0.194
  C_L 0.7 / 1.3               3.091 +- 0.201 / 2.943 +- 0.309
  tau 2 s / 12 s              3.023 +- 0.239 / 2.996 +- 0.273
  => the width of the spin prior is the only knob that matters at the 0.1 m level;
     C_D and tau are irrelevant at that level.

## Observation-noise sweep (entry y, sigma_pos in mm)
  E1: 2.9 -> 2.495 +- 0.166 (P_in 0.992); 5 -> 0.287 (0.921); 10 -> 0.594 (0.726);
      15 -> 0.966 (0.580); 20 -> 1.393 (0.470); 30 -> 2.239 (0.333). mean y constant.
  E2: 5 -> 2.815 +- 0.189 (0.969); 10 -> 3.002 +- 0.266 (0.891);
      15 -> 3.104 +- 0.319 (0.805); 20 -> 3.155 +- 0.319 (0.807).
  Tracking noise measured from open data (local cubic fits, 9-frame windows, 25 Hz):
      IDSSE 0.0558 m (x), 0.0544 m (y), 16791 / 15096 windows
      Metrica 0.0651 m (x), 0.0640 m (y)
  Hawk-Eye / FIFA GLT acceptance: 1.5 cm position accuracy, goal signal within 1 s
  (FIFA GLT test manual 2020, s4.1.2.2-4.1.2.4).  EPTS ball tracking is certified in 2-D only.

## Design study (linearised, sigma_pos = 10 mm, priors as above)
  sigma_y (mm) for n fixes evenly spaced at 0.05, 0.10, ..., n*0.05 s:
    n=2: 288 (sigma_pos 1 mm: 50; 5 mm: 202; 20 mm: 338; 50 mm: 388)
    n=3: 147 (19 / 81 / 228 / 304)
    n=4: 72  (10 / 37 / 130 / 223)
    n=6: 21  (3 / 11 / 40 / 92)
    n=8: 15
  If spin is known to 0.6 rad/s: n=2, sigma_pos 5/10/15/20/50 mm -> 18/35/51/68/170 mm
  Best two-fix timing (grid 0.03-0.31 s): t1 = 0.12 s, t2 = 0.31 s -> sigma_y = 22 mm
  (13x better than 0.05/0.10 s at the same 10 mm noise). The information is in the *baseline*.
  Rule of thumb: sigma_y ~ sigma_pos * 11 / (v Dt) with Dt the fix separation -> the
  dominant lever on the answer is the separation of the two observations, not their precision.

## Validation (synthetic, 120 trials per truth, sigma_pos = 10 mm) -- FILL FROM results/valid.json
  see results/valid.json; E1 bias/RMSE/coverage vs E2 bias/RMSE/coverage for physics-truth
  and kinematic-truth ensembles.

## Timing / feasibility (for the discussion)
  A MAP + Laplace-covariance prediction takes 8 ms + 0.5 ms on one core (MCMC 2.7 s offline);
  the two fixes arrive at 0.10 s, so the point estimate and its covariance are available
  by ~0.11 s after the kick, i.e. 0.28 s before the ball reaches the plane (0.39 s) and
  0.16 s before a keeper must commit (115 ms before contact).
