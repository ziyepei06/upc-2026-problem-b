# Main-line specification (authoritative for the rewrite)

The paper is re-based on the DS/GPT analysis line. `FINAL_NUMBERS.md` in the workspace
repo root (`../FINAL_NUMBERS.md` in the working tree) is the single source of
truth for every number. This file adds the writing rules, the figure/table inventory and
the page budget. If a number appears in both this file and `FINAL_NUMBERS.md`, they agree;
if `FINAL_NUMBERS.md` has more detail, use it.

## 1. Headline statements (locked)

* Point estimate of the crossing: **(y, z) = (2.7, 1.95) m**, reached at
  **t = 0.399 s** after contact.
* Probability that the whole ball passes inside the frame: **P(goal) ≈ 0.97**
  (best estimate 0.972 ± 0.0004; estimator range 0.965–0.997; scenario range 0.87–1.00).
* Main uncertainty statement: a **model envelope** `y ∈ [2.14, 3.24] m`,
  `z ∈ [1.62, 2.27] m`. It is a range of plausible predictions under different model
  forms and priors, **not** a credible interval — say so in the same sentence.
* Bayesian conditional estimate (the 6-parameter fitted posterior, GPT v3):
  y = 2.7847 m (90 % [2.45, 3.12]; 95 % [2.392, 3.188]),
  z = 1.9452 m (90 % [1.62, 2.27]; 95 % [1.555, 2.340]),
  t = 0.3988 s, P(goal) = 0.97225 ± 0.00043, N = 131 072 draws, ESS = 109 585.

## 2. Model ladder (Table 2 of the paper)

| Model | k | Physics | rms (mm) | y (m) | z (m) |
|---|---|---|---|---|---|
| M0 constant acceleration | 6 | one parabola per coordinate, exact | 0.000 | 2.4948 | 1.9204 |
| M1 gravity only | 3 | a = −g ẑ | 10.376 | 3.2170 | 1.9570 |
| M2 gravity + drag | 4 | quadratic drag, no lateral force | 5.029 | 3.2172 | 1.8790 |
| M3 drag + lateral | 5 | drag plus constant lateral acceleration | 0.542 | 2.5486 | 1.8931 |
| M4 force decomposition | 0 | drag + Magnus read off M0 (derived) | — | 2.6481 | 1.9568 |
| M5 joint Bayes, MAP | 7 | drag + Magnus, prior-weighted | 0.869 | 2.7073 | 1.9541 |
| M6 transverse fit | 6 | drag + two transverse force coefficients | — | 2.7847 | 1.9452 |

Supporting numbers: M0 gives **a = (−16.000, −16.000, −12.000) m s⁻²**,
**v0 = (30.4, 9.4, 6.9) m s⁻¹**, t = 0.405008 s. M4 gives C_D = 0.378,
C_L = 0.211, |a_M| = 11.05 m s⁻², ω̂ = (0.2357, −0.1035, −0.9663), t = 0.4007 s.
M3 gives v0 = (30.405, 9.410, 6.960) m s⁻¹, C_D = 0.319, a_y = −11.422 m s⁻².
M6 gives C_D = 0.3466, C_L = 0.159. M5 is quoted only as a prior-sensitivity aside.

## 3. Intervals, calibration, identifiability

* M3 profile 90 %: y [2.144, 2.954], z [1.860, 1.926]; Wald 90 %: y [2.108, 2.989],
  z [1.851, 1.935]. Curvature check σ_prof/σ_Wald = 1.001 (y) and 1.006 (z); M2's 0.572
  is the √3 of a 2-degree-of-freedom scale estimate, not non-quadraticity.
* Calibration with n = 1000 synthetic truths: coverage y 0.917, z 0.446;
  λ_y = 0.9605, λ_z = 4.9402; calibrated M3 90 % y [2.125, 2.972], z [1.685, 2.101].
  n = 50 single seeds scatter 3.9–5.8 in λ_z — never quote 3.88.
* Identifiability of the 7-parameter model: data Jacobian 6×7, rank 6; singular values
  327.227, 318.306, 217.813, 1.766, 1.765, 1.750; prior shrinkages of the transverse
  coefficients 8.69 % (along the flight direction), 69.99 %, 72.65 %; effective number of
  data-determined parameters ≈ 5.27. This is why the paper fits six parameters.
* Validation with n = 500: correct-model truths cover y 0.918 [0.891, 0.939],
  z 0.928 [0.902, 0.948], 2-D 0.852 [0.818, 0.880], mean 90 % half-widths
  (0.497, 0.0477) m; deliberately wrong (vertical-Magnus) truths cover y 0.872,
  z 0.338, 2-D 0.292, half-widths (0.523, 0.0501). Note the 0.8 m s⁻¹ launch-velocity
  jitter of this design.

## 4. Assumption-scenario envelope (main-line robustness, 26 single-factor scenarios)

* P(goal) ranges 0.735 (σ = 20 mm) to 0.9997 (σ = 2.5 mm).
* Scenario means: y 2.687–3.130 m, z 1.846–2.042 m.
* Union of the marginal 95 % intervals: y [2.259, 3.935] m, z [1.169, 2.698] m.
* Four joint-stress scenarios: P(goal) 0.561–0.987, y mean 2.733–3.339 m.
* Timing jitter, correlated noise, wind, drag-level and prior-width variants are all
  inside this set. Wind and launch velocity are not separable from two fixes: after a
  refit a 3 m s⁻¹ crosswind moves the entry point by **0.007 m**, not 0.11 m.
* Explicitly: this is a discrete scenario envelope, not a confidence or credible interval.

## 5. Control models from our own earlier analysis (comparison only)

* E1, constant acceleration from the same fixes: y = 2.50 ± 0.60 m,
  z = 1.93 ± 0.60 m, P(goal) = 0.73.
* E2, our drag-plus-Magnus Bayesian fit with literature priors (C_D = 0.197 ± 0.029,
  |ω| = 52 ± 22 rad s⁻¹): y = 3.03 ± 0.28 m, P(goal) = 0.89.
* Equal-weight average: y = 2.76 ± 0.54 m, P(goal) = 0.81.
* These are reported once, as a comparison, and they are what motivated the model
  envelope; they are not the paper's answer.

## 6. Forbidden numbers and phrases

Do not use: "100 repetitions", 3.88, 0.970 as an exact value, [2.2, 3.05], [1.4, 2.25],
"conservative rather than biased", "0.11 m" for the wind, or any of our mixture numbers
as a headline. σ_pos = 5 mm must be described as half of the last printed digit of the
given coordinates, not as a camera calibration.

## 7. Figures (files in `figures/`, all PDF)

| label | file | content |
|---|---|---|
| fig:geometry | f01-geometry.pdf | geometry, both fixes, goal plane |
| fig:trajectories | f02-trajectories.pdf | M2, M3 and M6 trajectories to the plane |
| fig:goalplane | f03-goalplane.pdf | 131 k posterior draws in the (y, z) plane |
| fig:marginals | f04-histograms.pdf | marginal densities of y and z |
| fig:identifiability | f05-identifiability.pdf | spectrum of the data Jacobian, ridge |
| fig:scenarios | f08-sensitivity.pdf | 26-scenario envelope |
| fig:noise | f09-noise.pdf | answer and P(goal) vs assumed noise |
| fig:design | f10-design.pdf | σ_y vs number of fixes and precision |
| fig:timing | f11-timing.pdf | two-fix timing map |
| fig:ladder | f13-models.pdf | model ladder M0–M6 |
| fig:aero | f14-aero.pdf | C_D(Re), C_L(S_p) |
| fig:timeline | f15-timeline.pdf | prediction window vs keeper deadline |

## 8. Page budget

Target 30–32 pages total: body §1–§9 about 22 pages, AI report ≤3 pages, references 1
page, appendix ≤4 pages. No ToC. Keep twelve figures and eight tables. Write compactly:
no bullet-point lists longer than four items, no repetition between §6 and §7.

## 9. Style

Formal academic English, third person or "we"; no Chinese characters anywhere; no
colloquialisms ("guesses", "a coin flip", "scouting report"); define every symbol at first
use; every number carries units and, where meaningful, an uncertainty; state the scope of
each claim in the same sentence as the claim.
