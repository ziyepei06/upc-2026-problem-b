# UPC 2026 Problem B — predicting the goal-plane entry of a penalty kick

**Team 230** · University Physics Competition 2026 · Problem B
[中文说明](README.zh-CN.md)

A penalty is given by the kick point $(0,0,0.11)$ m and two measured ball
positions, at $t=0.05$ s and $t=0.10$ s after the kick. This repository holds the
paper (LaTeX source and compiled PDF), the code that produced every number in
it, and the machine-readable results.

## Result

A six-parameter transverse-force model — launch velocity $\vect v_0$ (3),
drag coefficient $C_D$ (1), and two coefficients of a force transverse to the
flight (2) — is fitted to the six coordinates by our own damped Gauss–Newton
routine at an assumed 5 mm position noise, and the posterior is propagated to
the plane $x=11$ m:

| quantity | value |
|---|---|
| entry point | $y=2.785$ m, $z=1.945$ m at $t=0.399$ s |
| position in the goal | 0.87 m inside the right post, 0.50 m below the crossbar |
| 90 % intervals | $y\in[2.455,3.121]$ m, $z\in[1.617,2.276]$ m |
| probability the whole ball passes cleanly | 0.971 ± 0.001 |
| model envelope (two complete models) | $y\in[2.125,3.121]$ m |
| scenario envelope (25 one-factor + 4 joint-stress cases) | $y\in[2.259,3.935]$ m |

What the paper adds beyond the number:

- **The second difference of the data is noise.** With two fixes 0.05 s apart the
  acceleration is uncertain by $\sigma_a=\sqrt6\,\sigma_{\rm pos}/\Delta t^2
  \approx 9$ m s$^{-2}$ at 5 mm, so the force that bends the shot is only weakly
  determined; the paper quantifies this with a Fisher analysis and a rank-6
  Jacobian whose three small singular values sit $10^{16}$ below the large ones.
- **Model form, not sampling, sets the width.** Two 90 % marginals are not a 90 %
  joint region (Appendix F), and the union over admissible models is five times
  the statistical interval.
- **Real tracking data cannot settle $C_D$.** Fitting the drag-only solution to
  350 one-second windows of two public 2-D datasets returns $C_D=0.38$
  (interquartile 0.28–0.73) on raw and $0.11$ on smoothed provider data, while
  the same estimator recovers $0.200\to0.209$ on synthetic flights: the spread
  belongs to the tracking pipeline, not to the ball (§6.8).
- **The two fixes are not a typical free flight.** They demand a horizontal
  deceleration of 22.6 m s$^{-2}$ (2.3g) at 31 m s$^{-1}$, above the 95th
  percentile of real trajectories.
- **The method transfers as an observation-design statement.** The problem's two
  epochs give $\sigma_y=211$ mm; a first fix at 0.11 s and a second at 0.31 s
  give 12 mm, a factor 18.

## What is in the repository

| path | contents |
|---|---|
| `230B.pdf` | the paper (36 pages, 15 figures, 10 tables, 26 references) |
| `main.tex`, `sections/`, `figures/`, `refs.bib` | LaTeX source; `make` builds it |
| `code/` | every script that produced a number: fit (`fit6.py`, `fit6_fix2.py`), sampler, design study (`design6.py`), validation (`valid6b.py`), real-data drag (`real_cd.py`), figure scripts, `predict.py` (the appendix listing), `audit.py` |
| `results/` | JSON output of each script, plus `chains/` |
| `notes/` | the numbers digest used while writing |

## Build

```bash
make          # xelatex, then the submission checks
make pdf      # also copies main.pdf to 230B.pdf
```

Requires TeX Live (XeLaTeX, `booktabs`, `siunitx`, `natbib`), Python 3 with
`numpy`, `matplotlib`, `pymupdf`. The figure scripts read the reference PDFs and
datasets through `UPC_ROOT`; the two that need no external data are
`code/mkfig_a.py` and `code/mkfig_c.py`.

## Numbers and how they were checked

`code/audit.py` recomputes the quoted quantities from the raw result files.
Two independent integrators (our RK4 and a DOP853 reference) agree to
$2\times10^{-6}$ m, two independent samplers of the same posterior to 3 mm, and a
second, independently written implementation of the pipeline reproduces the
answer to 4 cm. The end-to-end coverage is measured on 750 synthetic truths:
500 on the five-parameter model, 50 on the adopted six-parameter model, and 25
under a law containing an extra vertical force.

## Data

The aerodynamic coefficients, the spin-decay range, the measured flight time and
the tracking precision are taken from the literature cited in the paper
(Carré et al. 2005; Goff & Carré 2012; Asai et al. 2020; Tsukada & Sakurai 2008;
Mehta 1985; Bassek et al. 2025; Metrica sample data). No public dataset contains
a penalty with a measured spin and a known goal-plane crossing, which is why the
end-to-end test is synthetic while the measurement half is real.

## Licence

Text, code and figures in this repository: MIT (see `LICENSE`). Third-party
papers and datasets remain under their own terms and are not redistributed here.
