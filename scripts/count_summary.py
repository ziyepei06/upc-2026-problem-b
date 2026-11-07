#!/usr/bin/env python3
"""Count the words in the Summary/Abstract and check the 300-word limit.

Usage:
    python3 scripts/count_summary.py paper.pdf            # count on page 1
    python3 scripts/count_summary.py main.tex --tex       # count before building
    python3 scripts/count_summary.py paper.pdf --limit 300 --verbose

The rule: page 1 holds the team number, the problem letter, the title and a
Summary of at most 300 words. Equations count as one word each; tables and
figures do not belong on page 1.
"""
import argparse
import os
import re
import sys

MAX = 300


def strip_latex(s: str) -> str:
    s = re.sub(r"(?m)(?<!\\)%.*$", "", s)                 # comments
    s = re.sub(r"\\begin\{(equation|align|gather)\*?\}.*?\\end\{\1\*?\}", " EQ ", s, flags=re.S)
    s = re.sub(r"\$[^$]*\$", " EQ ", s)                   # inline math
    s = re.sub(r"\\(cite|ref|eqref|label)\{[^}]*\}", " CIT ", s)
    s = re.sub(r"\\(textbf|textit|emph|texttt|mathrm|text|SI|si|num|qty)\{([^{}]*)\}", r"\2", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", s)    # remaining commands
    s = s.replace("{", " ").replace("}", " ").replace("~", " ")
    return s


def without_keywords(text: str) -> str:
    """Drop the keyword line: it sits under the abstract but is not part of
    the 300-word limit (the rule counts the summary; keywords are a label)."""
    keep, skipping = [], False
    for line in text.splitlines():
        if re.match(r"(?i)^\s*(?:\\noindent\s*)?(?:\\textbf\{)?key\s*words?\s*:", line):
            skipping = True
            continue
        if skipping:
            if not line.strip():
                skipping = False
            continue
        keep.append(line)
    return "\n".join(keep)


def words(text: str) -> int:
    return len([w for w in re.split(r"\s+", text.strip()) if re.search(r"[A-Za-z0-9]", w)])


def expand_inputs(path: str, seen: set | None = None, depth: int = 0) -> str:
    """Read a .tex file and inline its \\input{} files (multi-file projects)."""
    seen = seen or set()
    path = os.path.abspath(path)
    if path in seen or depth > 10:
        return ""
    seen.add(path)
    src = open(path, encoding="utf-8", errors="ignore").read()
    base = os.path.dirname(path)

    def repl(m):
        name = m.group(1).strip()
        if not name.endswith(".tex"):
            name += ".tex"
        for cand in (os.path.join(base, name), os.path.join(base, "sections", name)):
            if os.path.exists(cand):
                return expand_inputs(cand, seen, depth + 1)
        return ""

    return re.sub(r"\\(?:input|include)\{([^}]*)\}", repl, src)


def from_tex(path: str) -> str:
    src = expand_inputs(path)
    # Two accepted layouts: the plain \section*{Summary} skeleton, and the
    # summary sheet in the multi-file template, where the heading is a centred
    # \scshape line inside the titlepage environment. "Abstract" is accepted
    # too: the rule calls the item a summary, but 64 of the 74 gold papers
    # label it Abstract.
    pats = [
        (r"\\section\*?\{\s*(?:Summary|Abstract)\s*\}(.*?)"
         r"(?=\\section|\\newpage|\\clearpage|\\end\{document\})", None),
        (r"\\(?:large|Large|LARGE|normalsize)\s*\\scshape\s*(?:Summary|Abstract)\s*\}"
         r"(.*?)(?=\\end\{titlepage\})", None),
    ]
    for pat, _ in pats:
        m = re.search(pat, src, flags=re.S)
        if m:
            body = m.group(1)
            # the summary sheet sets its measure right after the heading; those
            # arguments are not prose and must not be counted
            body = re.sub(r"\\(?:setlength|setstretch|addtolength)\{[^{}]*\}\{[^{}]*\}", " ", body)
            body = re.sub(r"\\[vh]space\*?\{[^{}]*\}", " ", body)
            body = re.sub(r"\\(?:setstretch|fontsize|vskip)\s*\{[^{}]*\}", " ", body)
            return strip_latex(without_keywords(body))
    sys.exit("could not find a Summary/Abstract block in " + path)


def from_pdf(path: str) -> str:
    try:
        import pymupdf
    except ImportError:
        sys.exit("need PyMuPDF: python3 -m pip install --user pymupdf")
    doc = pymupdf.open(path)
    page = doc[0]
    txt = page.get_text()
    m = re.search(r"(?im)^\s*(?:summary|abstract)\s*$", txt)
    body = txt[m.end():] if m else txt
    # drop a running header on page 1 (team number / problem letter)
    lines = [ln for ln in body.splitlines()
             if not re.match(r"(?i)^\s*(team\s*\d+|problem\s*[ab])\s*$", ln)
             and not re.match(r"^\s*\d{1,3}\s*$", ln)]
    return without_keywords("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("--tex", action="store_true", help="input is a .tex source")
    ap.add_argument("--limit", type=int, default=MAX)
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    if not os.path.exists(a.input):
        sys.exit("no such file: " + a.input)

    raw = from_tex(a.input) if a.tex else from_pdf(a.input)
    n = words(raw)
    status = "OK" if n <= a.limit else "OVER by %d" % (n - a.limit)
    print("summary words: %d / %d  [%s]" % (n, a.limit, status))
    if a.verbose:
        print("-" * 60)
        print(raw.strip())
        print("-" * 60)
        print("structure check:")
        low = raw.lower()
        for label, pats in [
            ("situation/approach", ("we ", "model", "approach", "method")),
            ("assumptions", ("assum", "neglect", "ideal", "ignore")),
            ("numbers with units", (r"\d",)),
            ("validation", ("valid", "compar", "agree", "experiment", "data")),
            ("uncertainty", ("uncertain", "error", "±", "+-", "range")),
        ]:
            hit = any(re.search(p, low) if p != r"\d" else re.search(p, raw) for p in pats)
            print("  [%s] %s" % ("x" if hit else " ", label))
    return 0 if n <= a.limit else 1


if __name__ == "__main__":
    sys.exit(main())
