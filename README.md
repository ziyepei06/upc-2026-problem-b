# UPC 2026 Problem B — predicting the goal-plane entry of a penalty kick

**Team 230** · University Physics Competition 2026 · Problem B
[中文说明](README.zh-CN.md)

A penalty kick is given by two measured ball positions — at $t=0.05$ s and
$t=0.10$ s after the kick, plus the kick point itself. This repository holds the
paper (LaTeX source and compiled PDF), the analysis code that produced every
number in it, and the machine-readable results.

## Result

Six coordinates determine the six acceleration components. We fit a
six-parameter transverse-force model — launch velocity $v_0$ (3), drag
coefficient $C_D$ (1) and two transverse coefficients (2) — with our own damped
Gauss–Newton routine at an assumed 5 mm measurement noise, and propagate the
posterior to the goal plane:

| quantity | value |
|---|---|
| entry point at $x=11$ m | $y=2.785$ m, $z=1.945$ m, $t=0.399$ s |
| position inside the frame | 0.87 m inside the right post, 0.50 m below the crossbar |
| 90 % intervals | $y\in[2.455,3.121]$ m, $z\in[1.617,2.276]$ m |
| probability the whole ball passes | 0.971 |
| model-form envelope (four complete models) | $y\in[2.125,3.121]$ m |

Main findings:

- **$C_D$ is not identifiable** from two fixes. Over $C_D=0.15$–$0.50$ the
  residuals stay flat (12–18 mm) and the entry point moves only in $z$ (≈5 cm);
  $y$ does not move at all.
- **Spin and lateral velocity are degenerate.** At the measured lift scale
  $C_L\approx0.7\,S_p$ a spin of 10 rad/s shifts the entry point by ≈7 cm;
  the data prefer a mild negative spin only at $\Delta\chi^2\approx1.4$
  ($p\approx0.24$) — not significant.
- **Model form dominates the answer.** The statistical interval is ≈3.5 cm; the
  envelope over physically admissible models is ≈1 m. The honest answer is a
  band, not a point.
- The pipeline is validated by injection–recovery on synthetic truths, and the
  paper includes a design study of what one additional tracking frame would buy.

## Repository layout

```
main.tex  preamble.tex  config.tex  refs.bib   LaTeX source; config.tex = team number + title
sections/                                       one file per section (00-frontmatter … 11-ai-report)
figures/                                        vector PDF figures used by the paper
code/                                           analysis code (numpy + matplotlib only)
results/                                        JSON outputs behind every number in the paper
notes/                                          derivations and the number-provenance index
scripts/                                        build/length/compliance helpers
230B.pdf                                        submission copy (team 230, problem B)
```

## Build the paper

```bash
make        # latexmk + the submission checks (page count, summary length, naming)
make pdf    # also copies the result to 230B.pdf
```

Needs a TeX Live with `latexmk`, `xelatex`/`pdflatex` and `bibtex`.

## Reproduce the analysis

```bash
cd code
python3 fit6.py      # damped Gauss–Newton fit + intervals   -> ../results/fit6.json
python3 an_valid.py  # synthetic injection–recovery
python3 mkfig_a.py   # figures that need no external data
```

The code needs only `numpy` and `matplotlib`: the fit, the importance sampler and
the figures were written without SciPy. Every number quoted in the paper is
mapped to its producing file in `notes/numbers.md`.

A few figure scripts also read the reference tracking datasets and the results
of the parallel pipeline used for the cross-check; those archives are not
redistributed here. Point `UPC_ROOT` at a directory that contains them:

```bash
export UPC_ROOT=$HOME/upc-2026-data
# <UPC_ROOT>/2026/UPC 2026 Problem B 资料库/08_数据集/...
# <UPC_ROOT>/GPT/prediction_revision_v5/results/...
```

## Data

The two measured positions and the geometry are from the problem statement. All
derived quantities are in `results/` (JSON). The external tracking datasets used
to calibrate the measurement-noise model (SoccerNet-v3D, SkillCorner sample,
StatsBomb open data) are cited in `notes/` and not redistributed.

## Use of AI

`config.tex` sets `\aiusedtrue`, and the section *Report on the Use of AI*
documents the tools used (DeepSeek Harness, model `deepseek-flash`) and the
division of labour between them and the team. Page 1 follows the UPC rule:
team number, problem letter, title and summary only — no names, no school.

## License and citation

MIT — see `LICENSE` and `CITATION.cff`.
