# notes/ — scratch space

Anything that is not part of the submitted PDF lives here. Nothing in this
folder is compiled, so there is no penalty for rough work.

Suggested files, all optional:

| File | Purpose |
|---|---|
| `numbers.md` | every number that appears twice (summary, results, conclusion), written once with its uncertainty |
| `derivation-<topic>.md` | algebra that is too long for the body but might move into the appendix |
| `data-<source>.md` | provenance of each data point: URL, table number, figure digitised or quoted |
| `params.csv` | nominal parameter values and their uncertainties, the source of Table 3 |
| `meeting-log.md` | decisions and who owns which section, one line per decision |

## The one rule

A number written in two places must be copied from one file — usually
`numbers.md` — and never retyped from memory. Divergent digits between the
summary and the conclusion are the most common self-inflicted deduction in
this contest, and they are invisible to everyone except the judge who
checks.

## Cleaning up

Nothing here ships with `000A.pdf`. If you keep the folder in a shared
repository, add `notes/` to the repository's `.gitignore` only if it
contains nothing you need to defend later — the meeting log usually is
worth keeping.
