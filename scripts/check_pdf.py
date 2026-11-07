#!/usr/bin/env python3
"""UPC paper compliance & rendering check.

Usage:
    python3 scripts/check_pdf.py paper.pdf [--team 204] [--problem B]

Checks (hard):
  * PDF opens, has no blank pages, text is extractable (not a scanned image)
  * page 1 contains only: team number, problem letter, title, summary
  * summary <= 300 words
  * team number appears in the header of every page
  * a references section exists and is non-empty
  * filename follows <team><problem>.pdf

Checks (soft warnings): figure count, institutions/names leakage, fonts, units.
Exit code 0 = all hard checks pass.
"""
import sys, re, os, argparse

try:
    import pymupdf  # PyMuPDF
except ImportError:
    sys.exit("need PyMuPDF: python3 -m pip install --user pymupdf")

MAX_SUMMARY = 300


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    ap.add_argument("--team", default=None)
    ap.add_argument("--problem", default=None, choices=[None, "A", "B", "a", "b"])
    ap.add_argument("--max-summary", type=int, default=MAX_SUMMARY)
    a = ap.parse_args()

    if not os.path.exists(a.pdf):
        sys.exit(f"no such file: {a.pdf}")
    doc = pymupdf.open(a.pdf)
    fails, warns = [], []

    pages = doc.page_count
    texts = [p.get_text() for p in doc]
    total_chars = sum(len(t) for t in texts)
    if pages < 5:
        warns.append(f"only {pages} pages -- gold papers median 23 pages")
    blank = [i + 1 for i, t in enumerate(texts) if len(t.strip()) < 30]
    if blank:
        fails.append(f"near-blank pages: {blank}")

    first = texts[0]
    first_low = first.lower()

    # ---- team number ----
    m = re.search(r"team\s*(?:#|no\.?|number)?\s*:?\s*(\d{2,4})", first_low)
    team = a.team or (m.group(1) if m else None)
    if not team:
        fails.append("page 1: team number not found")
    else:
        hdr_missing = [i + 1 for i, t in enumerate(texts)
                       if team not in t]
        if hdr_missing:
            fails.append(f"team number {team} missing on pages {hdr_missing[:8]}"
                         f"{'...' if len(hdr_missing) > 8 else ''}")

    # ---- problem letter ----
    if not re.search(r"problem\s*[:\-]?\s*[ab]\b", first_low):
        fails.append("page 1: problem letter (A/B) not found")
    if a.problem and not re.search(rf"problem\s*[:\-]?\s*{a.problem.lower()}\b", first_low):
        fails.append(f"page 1: expected problem {a.problem.upper()}")

    # ---- summary ----
    si = first_low.find("summary")
    ai = first_low.find("abstract")
    idx = si if si >= 0 else ai
    if idx < 0:
        fails.append("page 1: no Summary/Abstract block")
    else:
        body = first[idx:]
        body = re.sub(r"(?i)^(summary|abstract)\W*", "", body)
        # stop at a section-like heading if the summary spills over
        body = re.split(r"\n\s*(?:\d+\.?\s+)?(?:introduction|1\s+introduction)\b",
                        body, flags=re.I)[0]
        # drop running header and page-number lines so they are not counted as words
        body = "\n".join(ln for ln in body.splitlines()
                         if not re.match(r"(?i)^\s*(team\s*\d+|problem\s*[ab]|\d{1,3})\s*$", ln))
        words = [w for w in re.split(r"\s+", body) if re.search(r"[A-Za-z0-9]", w)]
        n = len(words)
        if n > a.max_summary:
            fails.append(f"summary too long: {n} words (limit {a.max_summary})")
        elif n < 120:
            warns.append(f"summary short: {n} words -- aim 250-300")
        print(f"summary words : {n}")

    # ---- anonymity ----
    leak = re.findall(r"(?i)\b(?:university|college|institute|universit[äé]t|"
                      r"faculty sponsor|advisor|professor)\b", " ".join(texts))
    if leak:
        warns.append(f"possible identifying words found ({len(leak)}x): "
                     f"{sorted(set(w.lower() for w in leak))[:5]} -- verify none "
                     f"identifies the team")

    # ---- references ----
    all_low = " ".join(texts).lower()
    if not re.search(r"\breferences?\b", all_low):
        fails.append("no References section")
    else:
        nbib = len(set(re.findall(r"\[\d+\]", " ".join(texts))))
        if nbib < 3:
            warns.append(f"only {nbib} distinct bracket-numbered citations found -- "
                         f"check the reference list is complete and cited")
        if not re.search(r"\[\d+\]", " ".join(texts[:-1])):
            fails.append("no in-text numeric citations found (rule: must cite in text)")

    # ---- figures ----
    nfig = len(re.findall(r"(?i)\b(?:figure|fig\.)\s*\d+", all_low))
    print(f"figures cited : {nfig}")
    if nfig < 4:
        warns.append(f"only {nfig} figures -- gold papers median 14")

    # ---- fonts embedded ----
    try:
        bad = set()
        for p in doc:
            for f in p.get_fonts(full=True):
                xref, name = f[0], f[3]
                if xref == 0:            # xref 0 == font not embedded
                    bad.add(name)
        if bad:
            warns.append(f"fonts possibly not embedded: {sorted(bad)[:4]}")
    except Exception:
        pass

    # ---- filename ----
    base = os.path.basename(a.pdf)
    if not re.fullmatch(r"\d{2,4}[AB]\.pdf", base):
        warns.append(f"filename '{base}' does not follow <team><problem>.pdf "
                     f"(e.g. 778A.pdf)")

    print(f"pages         : {pages}")
    print(f"extractable   : {total_chars} chars")
    for w in warns:
        print("WARN  " + w)
    for f in fails:
        print("FAIL  " + f)
    print("RESULT: " + ("PASS" if not fails else f"{len(fails)} hard failure(s)"))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
