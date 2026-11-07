# Evidence brief — 2026 UPC Problem B (penalty trajectory prediction)

Scope: every number below was read out of the resource library
`2026/UPC 2026 Problem B 资料库/`. Citation form: `folder/file.pdf p.N` = **PDF page N**
(printed journal page given separately where the offset matters). Nothing from
`problemB/`, `problemB_v2/`, `lit/`, `research/` was opened. Numbers marked
**[computed]** were calculated here from source numbers or from the problem data;
they are not quotations.

Problem data (for reference): origin at penalty mark, goal plane `x = 11 m`,
`−3.66 < y < 3.66`, `0 < z < 2.44`; ball centre at `t=0` is `(0, 0, 0.11)`;
`t=0.05 s → (1.50, 0.45, 0.44)`; `t=0.10 s → (2.96, 0.86, 0.74)`.
`z0 = 0.11 m` implies the problem assumes **ball radius 0.11 m, diameter 0.22 m**,
which matches the measured 218 mm ball of Goff & Carré and Carré's 218 mm full-size ball.

---

## A. Ball and pitch standards

Source: `11_官方规则与规格/IFAB-Laws-of-the-Game-2025-26.pdf`

| Quantity | Value | Page |
|---|---|---|
| Circumference | 68–70 cm (27–28 in) | p.26 §2 |
| Mass at start of match | 410–450 g | p.26 §2 |
| Pressure | 0.6–1.1 atm = 600–1100 g/cm² **at sea level** | p.26 §2 |
| Distance between posts (inside) | 7.32 m (8 yd) | p.24 §1 |
| Lower edge of crossbar | 2.44 m (8 ft) | p.24 §1 |
| Penalty mark from goal line (centre of mark to back edge of line) | 11 m (12 yd) | pp.22–23, §1 |
| Penalty area | 16.5 m deep, 40.3 m wide | p.22 figure/§1 |

Air density for match conditions: **NOT FOUND in IFAB**. Use the measured values
from the aerodynamics papers: ρ = 1.2 kg/m³ (`04_/GoffCarre2012-轨迹分析反演气动系数.pdf` p.2;
same value in `13_/2020-ISEA-旋转足球飞行轨迹与Magnus力实测.pdf` p.2) and
ρ = 1.23 kg/m³ (`04_/Nathan2008-棒球自旋马格努斯力.pdf` p.2).
ν = 1.54×10⁻⁵ m²/s (Goff & Carré p.2). m = 0.424 kg, A = 0.0375 m² (Goff & Carré p.2);
weight = 4.2 N (Goff & Carré p.2, "≈4.2 N"). Note IFAB's 410–450 g brackets 424 g.

---

## B. Drag: Cd(Re), drag crisis, surface effects

**Equations (all from `04_/GoffCarre2012-轨迹分析反演气动系数.pdf` p.2):**
`Re = vD/ν`, `F_D = ½ C_D ρ A v²`, `F_L = ½ C_L ρ A v²`, `S_p = Rω/v`.
Convenient scaling printed there: for D = 0.218 m, `Re×10⁻⁵ ≈ v/(7 m/s) ≈ v/(16 mph)`.

**Drag crisis (the single most important aerodynamic fact for us):**

| Claim | Value | Source, page |
|---|---|---|
| Cd range across transition | 0.5 → 0.2 | Carré2005 p.2 (abstract) |
| Crisis interval | Re ≈ 9×10⁴ – 1.3×10⁵ | Carré2005 pp.6–7 |
| Same interval in speed (full-size ball) | 6.1 m/s → 9.1 m/s | Carré2005 pp.6–7 |
| Crisis speed (free-kick trajectory fit) | ≈ 14 m/s | GoffCarre2012 p.5 |
| Cd at Re = 3.9×10⁵ (26 m/s central kick) | 0.22 | Carré2005 p.11 |
| Cd at Re = 2.7×10⁵ (18.5 m/s curled kick) | 0.22 | Carré2005 p.12 |
| Cd, non-spinning, Re = 2.23–2.76×10⁵ | ≈ 0.17 | GoffCarre2012 p.5 |
| Smooth sphere Cd, Re = 7×10⁴–5×10⁵ | 0.53 ± 0.02 (no crisis) | Carré2005 p.10 |
| Separation angle below/above crisis | Φ ≈ 198° (v=7.65 m/s) → Φ ≈ 97° (v=19.1 m/s) | GoffCarre2012 p.5 |
| Cd repeatability (force balance) | ±1.7 % on Cd (p.10); 0.01 N (0.5 %) on spinning drag (pp.9–10); 0.007 N (2.1 %) on Magnus (p.8) | Carré2005 |
| Force balance resolution | 0.018 N over 0–10 N | Carré2005 p.5 |

**[computed] (v, Re) for D = 0.22 m, ν = 1.54×10⁻⁵ m²/s:**
v = 6.1 → 8.7×10⁴; 9.1 → 1.30×10⁵; 15 → 2.14×10⁵; 20 → 2.86×10⁵;
25 → 3.57×10⁵; 30 → 4.29×10⁵; 35 → 5.00×10⁵.
So a penalty struck at 25–35 m/s sits **entirely above the drag crisis**, where
Cd is nearly flat. Quoted Cd in that band: 0.22 (Carré, 2.7–3.9×10⁵) and 0.17
(Goff & Carré, 2.2–2.8×10⁵). **[computed]** the problem data give mean speed
30.9 m/s → Re = 4.4×10⁵ → Cd ≈ 0.19 by log-log interpolation of the table in
`09_鹰眼系统_分析与模型/model/physics.py` (that file's Cd table alone is:
2×10⁴→0.50, 4×10⁴→0.48, 6×10⁴→0.44, 7×10⁴→0.40, 8×10⁴→0.35, 9×10⁴→0.30,
1×10⁵→0.25, 1.3×10⁵→0.21, 2×10⁵→0.20, 4×10⁵→0.19, 1×10⁶→0.18).

**The safe statement for the paper:** at 25–35 m/s, Cd = 0.17–0.22 for a modern
32-panel ball, with the crisis safely below (so small speed errors do not move Cd).
`09_.../notes/07-计算推理模型.md` line 56 uses Cd = 0.189 at 30.9 m/s and reports
that drag alone changes flight time 0.3253 → 0.3440 s and lowers impact height
1.873 → 1.836 m. `09_.../README_索引.md` and `notes/02` also assert 500 fps for the
GLT cameras (secondary — see warning in J).

**Surface / seam effects (Carré2005 p.2, pp.6–13):** seams trip the boundary layer
the same way dimples do, which is why a 32-panel ball has its crisis at
Re ≈ 1×10⁵ instead of the smooth-sphere value (no crisis up to 5×10⁵ in these
measurements). Panel roughness changes the *critical* Re, therefore changes Cd only
near the crisis; above 25 m/s the panel details matter mainly through the Magnus
coefficient. Mehta 1985 (`04_/Mehta1985-运动球类空气动力学综述.pdf` pp.29–30,
journal pp.179–180) gives the mechanism in words: post-critical spin delays
separation on the retreating side, advances it on the advancing side; dimples/seams
lower the critical Re. Mehta is an **OCR'd scan with garbled text and no
football-specific Cd table** — cite only for terminology/mechanism.

---

## C. Magnus / spin

**Force law:** `F_M = ½ C_L ρ A v²` directed along `ω × v` (GoffCarre2012 p.2;
`13_/2020-ISEA...pdf` p.2 writes the side force as `S = ½ C_s ρ U² A`).
Spin parameter `S_p = Rω/v` (GoffCarre2012 p.2; Carré2005 uses the same definition).

**Measured C_L (or C_s) vs S_p — do not use a single coefficient:**

| S_p | C_L / C_s | Implied C_L/S_p | Source, page |
|---|---|---|---|
| 0.43 (ω = 64 rad/s, v = 18.5 m/s) | C_m = 0.32 | 0.74 | Carré2005 p.12 |
| 0.50 (ω = 63.6 rad/s = 608 rpm, v ≈ 14 m/s) | C_L ≈ 0.35, "levels off" beyond | 0.70 | GoffCarre2012 p.5 |
| 0.18 (time-averaged) | 0.21 | 1.17 | Asai2020 p.6 |
| 0.22 | 0.24 | 1.09 | Asai2020 p.6 |
| 0.34 | 0.29 | 0.85 | Asai2020 p.6 |
| 0.14 / 0.18 / 0.23 (wind tunnel) | 0.20 / 0.23 / 0.24 | 1.43 / 1.28 / 1.04 | Asai2020 p.6 |

So **C_L/S_p ≈ 0.7–1.2** for a football in the range S_p = 0.15–0.5, falling as S_p
grows beyond ~0.3. `09_.../model/physics.py` uses the conservative linear
`C_L = 0.505 S_p` (attributed there to "Goff & Carré 2010"); **that exact 0.505
coefficient is not printed in any PDF in this library** — it is a modelling choice,
not a citable measurement. If we quote a single number, say "C_L ≈ 0.2–0.3 for
S_p ≈ 0.2–0.35 (Asai2020 p.6)".

**Spin rates that actually occur:**

| Situation | Spin | Source |
|---|---|---|
| Penalty-type kicks (30 kicks, 18 m, 22 analysed) | 2.3–10.1 rev/s (14–63 rad/s) | Tsukada2008 p.1 |
| Central power kick (Asai data, 26 m/s) | 1.9 rev/s (12 rad/s), S_p = 0.05 | Carré2005 p.11 |
| Curled/off-centre kick (18.5 m/s) | 10.2 rev/s (64 rad/s), S_p = 0.43 | Carré2005 p.12 |
| Measured curled kick (23 m to goal, 19 m/s) | 5.8 rev/s (36 rad/s) | Asai2020 p.2–3 |
| Max in wind-tunnel spin tests | ±2800 rpm (=46.7 rev/s); the paper prints "(150 rad/s)", which is an internal inconsistency — 2800 rpm is 293 rad/s. Do not quote 150 rad/s. | Carré2005 p.8 |
| Baseball for contrast | 1500–4500 rpm = 25–75 rev/s | Nathan2008 p.1 |

**Magnus force magnitude, real numbers:** Asai2020 p.3 — side force on a 19 m/s,
5.8 rev/s kick falls ≈ 2 N → ≈ 1 N over 1 s (−50 %), while C_s falls only 0.27 → 0.24
(−10 %); the fall is mostly the v² factor. For scale: ball weight = 4.2 N
(GoffCarre2012 p.2), so the Magnus force is **≈ 25–50 % of weight** — big enough to
move the impact point by ~0.7 m over 1 s (Asai2020 p.3: > 0.7 m horizontal
deflection in 1 s).

**Spin decay — the citation set for the decay term:**

| Claim | Value | Source |
|---|---|---|
| Spin loss over 18 m of flight | 5–15 % (linear speed loses 16–28 % over the same 18 m) | Tsukada2008 p.1 |
| Spin loss, 0.8 s | 5.8 → 5.7 rev/s (≈ 2 %) | Asai2020 p.5 |
| Spin loss during the modelled flight | "only 1 % in the worst case" | Carré2005 p.12 |
| Spin-axis direction | essentially unchanged over the flight | Tsukada2008 p.1 |
| Exponential time constant used in the library's own chain | τ = 6 s (reproduces 5–15 % per 18 m) | `09_.../model/physics.py` |
| Deceleration comparison | spin deceleration is smaller than linear deceleration and independent of initial spin | Tsukada2008 p.1 |

No paper prints a fitted τ; **[computed]** 5–15 % per 18 m at 25–30 m/s
(0.6–0.72 s flight) implies τ ≈ 4–12 s, so τ = 6 s is a defensible nominal with a
factor-2 uncertainty. **This is exactly the kind of parameter the sloppy-model
literature says we cannot pin down (see G).**

---

## D. Penalty kinematics in reality

| Quantity | Value | Source, page |
|---|---|---|
| Mean flight time, 88 penalties (elite takers) | 525 ms, SD 64 ms | `12_/2020-FrontPsychol...pdf` p.7 |
| Ball flight times quoted for penalties | "vary between 500 ms and 600 ms" | `12_/2025-arXiv-Stop-Guessing...pdf` p.2 |
| Flight time over 18 m, power kick / curled kick | 0.9 s / 1.4 s | Carré2005 p.12 |
| Ball speed used by the 2022 gold-medal paper B-358 | v0 = 25 m/s, ω = 15 k̂ rad/s | `01_/2022-GM-B-358.pdf` p.1 |
| Shot speeds in IDSSE sample (11–14.6 m from goal) | 9.6–19.6 m/s peak | `08_/IDSSE-J03WPY-shots.csv` |
| Height distribution, N = 638 penalties (low 0–81.3 cm, mid 81.4–162.7, high 162.8–244) | 31.3 % low, 36.0 % mid, 32.7 % high | `12_/2020-FrontPsychol...pdf` p.5, Table 1 |
| Save rate by height (33/200, 77/230, 19/208) | 16.5 % low, 33.5 % mid, 9.1 % high | same, p.7 |
| Older sample, 311 kicks | 56.6 / 30.4 / 12.9 % low/mid/high | quoted in same paper p.8 from Bar-Eli & Azar 2009 (**secondary — not in library**) |
| Mean lateral ball position at goal line | −227 cm (SD 64) left, +233 cm (SD 81) right; goal half-width 366 cm | same, p.7 |
| Mean ball height at goal line | 37 cm (SD 24) low kicks, 203 cm (SD 27) high kicks | same, p.7 |
| Top vs bottom outcome (empirical, large sample) | top shots missed 12.88 % vs bottom 2.48 %; top saved 5.30 % vs bottom 23.51 % | `12_/2025-Cloud...pdf` p.~7 |
| Equilibrium conversion | corner 81.4 %, centre 79.6 % (simulation, n = 10⁹) | same, p.6 |

Take-away for the paper: a real penalty crosses the goal plane at
**|y| ≈ 1.5–2.4 m, z ≈ 0.3–2.0 m**, i.e. the ball is *not* near the centre; and
**~1/3 of kicks are above 1.63 m**, so z-accuracy matters as much as y-accuracy.
The "top corner" target that the 2022 gold papers optimise is also where shots miss
the frame most often (12.88 % vs 2.48 %).

**[computed] data-consistency check on the problem's own three points** (this is a
finding, not a literature claim): v(0→0.05) = (30.0, 9.0, 6.6) m/s, |v| = 32.01 m/s;
v(0.05→0.10) = (29.2, 8.2, 6.0) m/s, |v| = 30.92 m/s; hence the implied acceleration
is **exactly (−16.0, −16.0, −12.0) m/s²** (the position increments are exact to the
printed 0.01 m). Compare drag + gravity at 32 m/s with Cd = 0.19: drag decel
10.5 m/s² along −v̂ plus g = 9.81 m/s² downward, i.e. about (−9.8, −2.9, −12.4) m/s².
The **z-component matches to 0.4 m/s², but x and y decelerate 3–5× more than drag
can explain**; a Magnus force is perpendicular to v and cannot do it either.
Interpretation for the team: either the test points come from a deliberately simple
constant-acceleration generator, or they contain a systematic offset. Either way the
two-frame vs three-frame choice is *the* modelling decision (compare
`09_.../notes/07` line 113: two-frame (3.118, 1.873) vs three-frame constant-accel
(2.617, 2.013), Δy = −0.50 m).

---

## E. Goalkeeper: the latency budget

| Quantity | Value | Source, page |
|---|---|---|
| Perceptual-motor delay (visuomotor) | 55–130 ms (Lee et al. 1983); 80–144 ms (Franks & Harvey 1997) | quoted in `12_/2020-FrontPsychol...pdf` p.10 |
| Horizontal dive component | already differentiated between left/right **115 ms before ball contact** (i.e. keepers commit before the kick) | same, p.7 |
| Vertical component | hand-height difference appears **245 ms after ball contact** | same, p.7 + abstract p.1 |
| Consequence | keeper uses kicker kinematics for left/right, ball trajectory for height | same, p.1, pp.9–10 |
| Ball flight time to react to | 500–600 ms (525 ± 64 ms measured) | Stop-Guessing p.2; FrontPsychol p.7 |

**Requirement this sets on us:** any prediction device must output a goal-plane
estimate **before ≈ 0.5 s after the kick**; if a keeper is to use it, effectively
by ~0.3 s (flight 0.525 s minus dive time). Our problem only gives data to 0.10 s,
so we have 0.10 s of observation and ~0.42 s of extrapolation → an extrapolation
factor of ~4× in time, and the entry point is ~8–11 m downstream.

---

## F. Tracking / measurement accuracy — the uncertainty anchor

### F.1 Goal-line technology (the only *official* position-accuracy number in football)

| Quantity | Value | Source, page |
|---|---|---|
| GLT positional accuracy required | **±1.5 cm** (target-wall and slider tests) | `06_/FIFA-门线技术测试手册2020.pdf` §4.1.2.2 p.16, §4.1.2.3 p.17, Tables 3–4 |
| Alert latency required | "GOAL" on referee's watch **within one second**, results ≤ 1 s | same, §4.1.2.4 p.17, Tables 5a/5b |
| Independence requirement | indication "immediate and automatically confirmed within one second" | same, p.iii |
| Test shots | target wall at 120/160/200/240 mm behind the line; B-shots at 50 km/h, C01 35 km/h, C02 70 km/h; 12 shots per B-row, 5 per C | same, Tables 3a/3b p.16 |
| High-speed verification | verification camera ≥ 30 fps on the goal line; "images captured at 2000 fps" for part of the tests | same, pp.29–31 |
| Hawk-Eye system | 7 cameras per goal, roof-mounted; **2 of 7 cameras suffice** to fix the 3-D position; all 7 add accuracy + redundancy; "millimetre accurate" | `06_/鹰眼GLT-官方说明.pdf` p.2; `06_/鹰眼GLT-工作原理.pdf` p.1 |
| Frame rate | **not stated in either Hawk-Eye doc**; 09_ notes say 500 fps for GLT cameras (secondary) | — see J |

### F.2 EPTS (player/ball tracking certification)

| Quantity | Value | Source, page |
|---|---|---|
| Reference system | 40 VICON cameras (34 Valkyrie + 6 Vantage V16), 40 × 30 m area, ≥35 × 25 m capture volume | `06_/FIFA-电子表现跟踪测试手册2025.pdf` §2.3.1 p.6 |
| Reference sampling | **all VICON trials at 100 Hz** | same, p.6 |
| Reference accuracy | calibration error verified **< 1 mm** | same, p.6 |
| Provider data handling | resampled to 50 Hz by linear interpolation; velocity 2nd-order low-pass Butterworth **1 Hz**; ball data **5 Hz** Butterworth; VICON reference speed by 3-point central difference | same, pp.8–9, §3/§3.6 |
| What is certified for the ball | only **2-D position (X, Y) and speed** — "no official z truth" | same, §5.3.5 p.26 |
| Position metric | straight-line difference between provider and VICON X/Y | same, §4.1 p.19 |
| Rating | 5-level colour scale in speed bands 0–5 … +25 km/h; thresholds learned from >1,000,000 points, >30 systems, 3 test events; **no published mm threshold** | same, §4.3 p.19 |
| Bandwidth hint | football tracking providers operate at 10 Hz and 25 Hz in the certification tests | same, p.9 |

### F.3 What this means for finite differences over Δt = 0.05 s **[computed]**

For independent position noise σ:

- velocity from a first difference: `σ_v = √2 σ / Δt`
- acceleration from a second difference: `σ_a = √6 σ / Δt²`
- jerk from a third difference: `σ_j = √20 σ / Δt³`

| σ (m) | Δt (s) | σ_v (m/s) | σ_a (m/s²) | σ_j (m/s³) |
|---|---|---|---|---|
| 0.015 (FIFA GLT requirement) | 0.05 | 0.42 | 14.7 | 537 |
| 0.005 | 0.05 | 0.14 | 4.9 | 179 |
| 0.002 | 0.05 | 0.057 | 2.0 | 72 |
| 0.015 | 0.01 (100 Hz) | 2.12 | 367 | 6.7×10⁴ |

Reading: **with ±1.5 cm positions, the acceleration is unmeasurable from two 0.05 s
intervals** (14.7 m/s² noise vs 12 m/s² signal — S/N ≈ 1). Even at σ = 5 mm the
acceleration is only 2σ. This is the quantitative justification for the paper's
choice between (i) assume constant velocity + physical drag/Magnus model and
(ii) fit a constant acceleration to three points and accept a huge error bar.
Converting a position error into an impact-point error: over the remaining flight
time T, a velocity error σ_v maps to `σ_x ≈ σ_v·T` (≈ 0.42 × 0.42 ≈ 0.18 m at
σ = 1.5 cm); an acceleration error maps to `σ_x ≈ ½ σ_a T²` (≈ 0.5 × 14.7 × 0.18
≈ 1.3 m, i.e. useless). **State this in the paper — it is the strongest argument
for not trusting three-point accelerations.**

Independent anchor for the position noise we should assume: `09_.../notes/07`
line 18 reports that a two-camera triangulation with 12 m baseline, f = 2000 px
and 1 px detection noise gives 7.7 mm reconstruction error (1σ per axis ≈ 8/22/15 mm).
So σ = 1–20 mm is the right envelope, and the FIFA ±1.5 cm figure is the
conservative end. Neither Hawk-Eye nor FIFA publishes a per-frame positional σ for
a *moving* ball in flight — both numbers are goal-line / 2-D metrics.

---

## G. Uncertainty and identifiability — how to talk about it

1. **Sloppiness (the framework to borrow).** `07_/Gutenkunst2007-sloppy模型的可证伪预测.pdf`
   (Ann. N.Y. Acad. Sci. 1115:203–211, 2007) p.1–2: for multiparameter nonlinear
   models the Hessian/Fisher eigenvalues span **more than 10⁶** (Fig. 1B, models with
   48 parameters / 68 data points included); parameters are well determined only
   along "stiff" directions; **predictions can be tightly constrained even when
   parameters are wildly uncertain** (p.2, §Constraining Parameters); linear
   (Gaussian) uncertainty propagation is "dangerous" and Monte-Carlo is the honest
   route (abstract, p.1). This is the exact justification for reporting a prediction
   interval on the impact point rather than error bars on Cd, C_L, τ.
2. **Limits of parameter estimation.** `07_/White2017-sloppy系统的参数估计局限.pdf`
   (arXiv:1602.05135v3): many parameters are "practically unidentifiable"; more
   importantly, **systematic model discrepancy can make a model less predictive after
   "optimal" parameter estimation** (abstract, pp.1–2, §Discussion). Use this to say
   that the two-frame/three-frame spread is a *model-discrepancy* term, not a
   statistical error.
3. **Spin is not observable from cheap tracking.** `13_/2023-arXiv-乒乓灰箱EKF自旋不可辨识.pdf`
   (Achterhold et al., L4DC 2023): spin "is hard to infer from position measurements
   of the trajectory" (p.3); their physics-based EKF (state 11-dim:
   p, v, ω, a_d, a_m; ODE `v̇ = −k_d|v|²v + k_m(ω×v) + g`) predicts the 1 s-ahead
   position to ~0.1–0.4 m, and **initialising ω from the launcher (rather than
   estimating it from positions) drastically reduces the 1 s prediction error**
   (Fig. 5, p.9; the plotted worst case is 0.393 m). Return rate 29/30 = 97.7 %
   (p.1). Translation to football: do **not** try to infer ω from three positions;
   either take ω from the kicker/impact model or marginalise over a prior and let
   the prior, not the data, carry the spin.
4. Practical recipe that follows: report `p_impact` with (a) a statistical band from
   position noise (F.3), (b) an envelope over the aerodynamic coefficient ranges
   (Cd 0.17–0.22, C_L/S_p 0.7–1.2, τ 4–12 s), (c) the two-frame/three-frame
   discrepancy as a separate systematic line. `09_.../notes/07` lines 108–123 does
   exactly this and is a legitimate in-library template.

---

## H. Datasets (all in `08_数据集/`) — verified by opening the CSVs

| File | Rows | Columns | Rate / units | 3-D? | Usable for? |
|---|---|---|---|---|---|
| `IDSSE-J03WPY-ball.csv` | 146,211 | `frame,time,x,y,speed,accel` | 40 ms steps = **25 Hz**; ISO-8601 time; x,y in **metres about pitch centre** (e.g. 6.40, −0.04); speed in km/h (e.g. 60.77); `accel` column is **all 0.00** | **No z** | continuous ball track; 25 Hz finite-difference noise study |
| `IDSSE-J03WPY-shots.csv` | 26 | `peak_speed_mps,distance_to_goal_m` | m/s and m | — | 26 shot summaries, 9.6–19.6 m/s at 11–14.6 m from goal (penalty-like range) |
| `IDSSE-J03WPY-events.xml` | — | event log (goals/shots with timestamps) | — | — | syncing shots to tracks |
| `ISSIA-3D.csv` | 2,465 | `'',frame,x_cam1..6,y_cam1..6,num_cameras,ball_3D,err_cam1..6,cam_list,opt_d_cam1..6,opt_e_cam1..6` | per-frame; `ball_3D` a string like `"[10.3869, -24.1457, -0.00766]"`; typically only **2 cameras** used (`cam_list "[3, 4]"`); reprojection errors ≈ 10–19 px | **Yes (x,y,z)** | genuine 3-D truth; also lets us propagate pixel error → 3-D error |
| `Metrica-第1场-球轨迹.csv` | 145,006 | `period,frame,time_s,ball_x,ball_y` | 0.04 s = 25 Hz; **normalised 0–1** pitch coordinates; NaNs present | No z | long continuous track, but normalised units must be rescaled with a known pitch size |
| `SoccerNet-v3D.csv` / `-parsed.csv` | 4,051 | `league,season,match,main_action,action,replay,img_w,img_h,is_ball,ball_bbox,calibration,JaC@0.005,JaC@0.01,JaC@0.02,ball_3D,rep_error,optimized_error,optimized_d,set` (+parsed x,y,z) | **one 3-D point per action/replay**, not a trajectory; `ball_3D` like `"[array([39.02, 2.3411, 0.0041256])]"`; rep_error 2.4–20.0, optimized_error 0.03–0.36 | Yes, single points | validation of single-frame 3-D localisation, NOT of trajectories |

Notes:
- **No dataset here contains a penalty kick with measured spin.** None has an ω
  column. Spin must come from the literature priors (§C) or from the kicker model.
- SoccerNet `main_action` is a numeric class code (values 0–23+), so "penalty"
  cannot be selected without the label mapping (not in the folder): treat as
  shots-in-general.
- The only dataset that can validate a *whole trajectory* in 3-D is **ISSIA-3D**
  (2 cameras per frame though — so it carries exactly the triangulation noise we are
  worried about). Metrica/IDSSE give long 2-D tracks at 25 Hz, ideal for testing the
  0.05 s finite-difference noise (§F.3) and drag fits (fit Cd from a 2-D track after
  projecting horizontally).
- `09_鹰眼系统_分析与模型/` (allowed: it is a system analysis, not a solution) contains
  `notes/01..08`, `model/{physics,estimator,kalman,triangulate}.py`, and a README
  index. Read `notes/02` for GLT acceptance numbers, `notes/04` for tennis spin
  bias, `notes/07` for the triangulation→KF→flight-model chain. Its physics.py
  holds the Cd table quoted in §B and the C_L = 0.505 S_p choice flagged in §C.

---

## I. Gold-medal papers: structure only (no physics copied)

### 2022 B (penalty-kick strategy) — `01_同题金奖论文_2022B/`

| Paper | Pages | Figures / tables | Refs | Skeleton |
|---|---|---|---|---|
| `2022-GM-B-358.pdf` | 24 | 9 figs / 1 table | 8 | Intro → Assumptions & Notation → Model (3.1 flight of spinning ball: buoyancy+gravity, drag, Magnus; 3.2 …) → numerical solution → reachable θ,φ ranges for upper corners → goalkeeper defence (vertical move / dive) → Monte-Carlo save probability → "dead angle" region → Appendix: full Python/`odeint` code |
| `2022-GM-B-471.pdf` | 44 | 10 figs / 9 tables | 11+ | Summary Sheet → Summary → Model Establishment & Analysis → Single-Objective Programming Model → Goalkeeper Save-Ball Model → Model Solving (**Euler method**) → Model Correction → Model Evaluation → Conclusion |
| `2022-GM-B-480.pdf` | 19 | 11 figs / 4 tables | 8 (list pp.13–14) | Intro → Restatement → Assumptions & Notations → Physical Analysis (4.1 ball motion, 4.2 initial velocities/spins for upper corner, 4.3 keeper defence, 4.4 best kick) → Model Evaluation (5.1 Strengths / 5.2 Weaknesses) → References → Appendix |

How each handled uncertainty + validation: **358** = Monte-Carlo over initial
conditions and keeper behaviour (p.~13 onwards), no empirical validation.
**471** = an explicit "Model Correction" section (error/robustness pass) and Euler
integration; validation is internal (parameter sweeps), not data.
**480** = purely analytic + strengths/weaknesses, no data validation.
→ Lesson: none of the three validates against real tracking data, and only 471 has a
dedicated correction section. That is a gap we can exploit (datasets in §H).

### 2025 B (external ballistics) — `02_同结构金奖论文_2025B_火炮弹道/`

| Paper | Pages | Figures / tables | Refs | Skeleton |
|---|---|---|---|---|
| `2025-GM-B-127.pdf` | 24 | 9 figs / 3 tables | 14 | Assumptions → Analysis → Notations → Modelling (simple case / air resistance / **uncertain target position** / **effects of wind** / numerics) → Results & Discussion (incl. horizontal distance variation, vertical height variation, wind) → Conclusion → **Strength & Weakness** → References → Appendix |
| `2025-GM-B-264.pdf` | 27 | 8 figs / 3 tables | 7 | … **Model Validation and Regression Fitting** → **Validation Methodology and Scheme** → **Validation Results Under Wind-Free Conditions** → **Quantitative Error Analysis** → conclusion |
| `2025-GM-B-188.pdf` | 69 | 8 figs / 3 tables | — | longest of the set; heavy appendix |

Full page/fig/table census of the folder (for calibration of length): 127→24 pp,
155→20, 188→69, 204→37, 207→38, 264→27, 311→17, 392→33, 428→32, 460→32 (text not
extractable — scanned), 496→27, 765→52, 796→111.

**Structural take-away for us:** the winning 2025 B papers separate
*model → validation against data → quantitative error analysis* into three named
sections (264), and treat *uncertain input* and *wind* as their own modelling
subsections (127). Copy that skeleton, not the physics: ours would be
Model → Two/Three-frame Identifiability → Uncertainty Propagation →
Validation on ISSIA-3D / IDSSE → Error budget vs FIFA ±1.5 cm.

---

## J. Warnings — scans, gaps, unverifiable claims

1. `04_/Mehta1985-运动球类空气动力学综述.pdf` is an **OCR'd scan**: text is garbled,
   figure captions truncated, no usable football Cd table. Cite only pp.29–30 for the
   separation mechanism and p.30 for spin-parameter/dimple statements.
2. `02_/2025-GM-B-460.pdf` — **text extraction fails on most pages (scanned images)**.
   Treat as unreadable; if we cite its structure, verify visually.
3. **Hawk-Eye frame rate (500 fps) is not in any primary document in the library.**
   It appears only in `09_.../notes/02` (secondary, sourced from EE Times) and in
   `09_.../README_索引.md`. The two official Hawk-Eye PDFs give only "7 cameras,
   millimetre accurate, 2-of-7 redundancy". If the paper needs a frame rate, quote
   the EPTS reference (100 Hz VICON) or GLT verification (2000 fps, ≥30 fps) instead,
   and attribute 500 fps as "vendor-reported, secondary".
4. **"< 1 mm" for Hawk-Eye** is from the 2026 vendor page as relayed by 09_ notes; the
   FIFA GLT manual's *requirement* is ±1.5 cm. Do not conflate them: ±1.5 cm is the
   acceptance criterion, <1 mm is a marketing/installation claim.
5. **C_L = 0.505 S_p** (used in `09_.../model/physics.py`) is not printed in any PDF
   here; the measured ratios are 0.7–1.2 (§C). If a referee asks, the honest citation
   is Asai2020 p.6.
6. **Magnus "curvature" numbers in Asai2020 p.3 are garbled in the PDF** (prints
   1.1 / 31.1 / 0.9 m⁻¹; the middle value is almost certainly 1.31). Do not quote
   31.1 m⁻¹.
7. **No library source gives air density for a specific match** (altitude/temperature/
   humidity). ρ = 1.2 kg/m³ is the papers' working value; a ±5 % density change is
   worth ~±2 % in drag force — mention it, don't model it.
8. Tsukada2008 is a **one-page conference abstract** (ISBS 2008, pp.698), 22 trials
   from 3 players — a weak but uniquely relevant source for spin decay. Do not
   overstate it as a full study.

---

## K. The 15 load-bearing numbers (with citation)

| # | Number | Citation |
|---|---|---|
| 1 | Cd drops 0.5 → 0.2 across the drag crisis | `04_/Carre2005-足球缝线与气动.pdf` p.2 |
| 2 | Crisis at Re = 9×10⁴–1.3×10⁵ = 6.1–9.1 m/s (full-size ball) | `04_/Carre2005...pdf` pp.6–7 |
| 3 | Cd = 0.22 at Re = 2.7×10⁵ and 3.9×10⁵ (modern ball) | `04_/Carre2005...pdf` pp.11–12 |
| 4 | Cd ≈ 0.17 at Re = 2.2–2.8×10⁵ (trajectory fit) | `04_/GoffCarre2012...pdf` p.5 |
| 5 | Hence Cd(30 m/s) ≈ 0.19 [computed, Re = 4.3×10⁵] — quote 0.17–0.22 | (3)+(4) |
| 6 | C_L ≈ 0.21/0.24/0.29 at S_p = 0.18/0.22/0.34 (⇒ 0.7–1.2 × S_p) | `13_/2020-ISEA...pdf` p.6 |
| 7 | Spin 5.8 rev/s at 19 m/s produces > 0.7 m sideways in 1 s; side force 2 N → 1 N | `13_/2020-ISEA...pdf` pp.2–3 |
| 8 | Spin decays 5–15 % per 18 m; axis unchanged; ⇒ τ ≈ 4–12 s, nominal 6 s | `04_/Tsukada2008...pdf` p.1; τ = 6 s in `09_.../model/physics.py` |
| 9 | Flight time of a real penalty = 525 ± 64 ms (N = 88) | `12_/2020-FrontPsychol...pdf` p.7 |
| 10 | Penalty speed used/measured: 25 m/s (gold-paper nominal), 19 m/s curled, 26 m/s power | `01_/2022-GM-B-358.pdf` p.1; Carré2005 p.12; p.11 |
| 11 | Where kicks cross the plane: 31.3 % low / 36.0 % mid / 32.7 % high; |y| ≈ 2.3 m lateral | `12_/2020-FrontPsychol...pdf` pp.5,7 |
| 12 | Goalkeeper visuomotor delay 55–144 ms; commits left/right 115 ms **before** contact, height 245 ms **after** | `12_/2020-FrontPsychol...pdf` pp.7,10 |
| 13 | FIFA GLT accuracy **±1.5 cm**, alert **≤ 1 s** | `06_/FIFA-门线技术测试手册2020.pdf` §4.1.2.2–4.1.2.4, pp.16–17 |
| 14 | EPTS reference = VICON **100 Hz**, calibration **< 1 mm**; ball certified only in 2-D; provider rates 10/25 Hz → 50 Hz | `06_/FIFA-电子表现跟踪测试手册2025.pdf` pp.6,8–9,26 |
| 15 | Finite-difference noise at Δt = 0.05 s: σ_v = 0.42 m/s, σ_a = 14.7 m/s² for σ = 15 mm [computed] — acceleration is unmeasurable | derived from (13) |
