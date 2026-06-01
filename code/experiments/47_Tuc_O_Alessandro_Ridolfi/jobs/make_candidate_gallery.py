#!/usr/bin/env python
"""
make_candidate_gallery.py

Generates a self-contained HTML gallery with all images embedded as base64.
The output file can be downloaded and opened in any browser on any machine.

Usage
-----
    python make_candidate_gallery.py \
        --shortlist-csv  <shortlist_coherent.csv> \
        --cand-dir       <stage3_viterbi/candidates/> \
        --fold-dir       <stage4_fold/> \
        --out-html       <candidate_gallery.html>
"""

import argparse
import base64
import csv
import glob
import os
import sys


def find_file(directory, pattern):
    matches = glob.glob(os.path.join(directory, pattern))
    return matches[0] if matches else None


def img_tag(abs_path):
    """Return an <img> tag with base64-encoded image, or a placeholder."""
    if abs_path and os.path.isfile(abs_path):
        with open(abs_path, "rb") as fh:
            b64 = base64.b64encode(fh.read()).decode("ascii")
        return '<img src="data:image/png;base64,{}" style="width:100%;border:1px solid #444;">'.format(b64)
    return ('<div style="width:100%;height:200px;background:#222;color:#888;'
            'display:flex;align-items:center;justify-content:center;'
            'font-size:12px;">Not found</div>')


def load_shortlist(csv_path):
    rows = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shortlist-csv", required=True)
    ap.add_argument("--cand-dir",      required=True)
    ap.add_argument("--fold-dir",      required=True)
    ap.add_argument("--out-html",      required=True)
    args = ap.parse_args()

    for p, label in [(args.shortlist_csv, "shortlist-csv"),
                     (args.cand_dir,      "cand-dir"),
                     (args.fold_dir,      "fold-dir")]:
        if not os.path.exists(p):
            sys.exit("ERROR: {} not found: {}".format(label, p))

    print("Loading shortlist: {}".format(args.shortlist_csv), flush=True)
    candidates = load_shortlist(args.shortlist_csv)
    print("  {} candidates.".format(len(candidates)), flush=True)

    def sort_key(r):
        try:
            return -float(r.get("snr_best", "nan"))
        except ValueError:
            return float("inf")
    candidates.sort(key=sort_key)

    n_total    = len(candidates)
    n_binary   = sum(1 for r in candidates if "binary"   in r.get("selection", ""))
    n_isolated = sum(1 for r in candidates if r.get("selection", "") == "isolated")
    n_known    = sum(1 for r in candidates if r.get("known_match", ""))

    print("Building HTML (embedding images as base64)...", flush=True)

    cards = []
    n_missing_spec = 0
    n_missing_poly = 0
    n_missing_kep  = 0

    for i, row in enumerate(candidates):
        cand_id      = row.get("cand_id",      "")
        freq_hz      = row.get("freq_hz",      "")
        snr_best     = row.get("snr_best",     "")
        path_std     = row.get("path_std_hz",  "")
        peak_ll      = row.get("peak_loglike", "")
        multiplicity = row.get("multiplicity", "")
        nt_values    = row.get("nt_values",    "")
        known_match  = row.get("known_match",  "")
        selection    = row.get("selection",    "")
        snr_poly     = row.get("snr_poly",     "")
        snr_kepler   = row.get("snr_kepler",   "")

        print("  [{}/{}] {} ...".format(i+1, n_total, cand_id), flush=True)

        is_known    = bool(known_match)
        card_border = "#c8a800" if is_known else "#2a2a2a"
        card_bg     = "#1a1500" if is_known else "#1a1a1a"

        spec_path = find_file(os.path.join(args.cand_dir, cand_id),
                              "*_spectrogram.png")
        poly_path = find_file(os.path.join(args.fold_dir, cand_id, "poly"),
                              "*.png")
        kep_path  = find_file(os.path.join(args.fold_dir, cand_id, "kepler"),
                              "*.png")

        if not spec_path: n_missing_spec += 1
        if not poly_path: n_missing_poly += 1
        if not kep_path:  n_missing_kep  += 1

        if "binary" in selection and "isolated" in selection:
            badge_col = "#4fc3f7"
        elif "binary" in selection:
            badge_col = "#81c995"
        else:
            badge_col = "#ce93d8"

        known_badge = ""
        if is_known:
            known_badge = (
                '<span style="background:#c8a800;color:#000;padding:2px 6px;'
                'border-radius:3px;font-size:11px;font-weight:bold;">'
                '{}</span>'.format(known_match)
            )

        sel_badge = (
            '<span style="background:{};color:#000;padding:2px 6px;'
            'border-radius:3px;font-size:11px;">{}</span>'.format(
                badge_col, selection)
        )

        try:
            freq_str = "{:.4f} Hz".format(float(freq_hz))
        except ValueError:
            freq_str = freq_hz

        try:
            snr_str = "{:.2f}".format(float(snr_best))
        except ValueError:
            snr_str = snr_best

        try:
            ll_str = "{:.1f}".format(float(peak_ll))
        except ValueError:
            ll_str = peak_ll

        try:
            std_str = "{:.4f} Hz".format(float(path_std))
        except ValueError:
            std_str = path_std

        card = """
<div id="{cand_id}" style="border:2px solid {card_border};background:{card_bg};
     border-radius:6px;padding:12px;margin-bottom:20px;">

  <div style="display:flex;align-items:center;gap:10px;margin-bottom:8px;flex-wrap:wrap;">
    <span style="font-size:16px;font-weight:bold;color:#eee;">#{rank} {cand_id}</span>
    {sel_badge}
    {known_badge}
    <span style="color:#aaa;font-size:13px;margin-left:auto;">
      <b style="color:#eee;">{freq_str}</b>
    </span>
  </div>

  <table style="font-size:12px;color:#ccc;border-collapse:collapse;
                margin-bottom:10px;width:100%;">
    <tr>
      <td style="padding:2px 10px 2px 0;"><b>S/N best:</b> {snr_str}</td>
      <td style="padding:2px 10px 2px 0;"><b>S/N poly:</b> {snr_poly}</td>
      <td style="padding:2px 10px 2px 0;"><b>S/N kepler:</b> {snr_kepler}</td>
      <td style="padding:2px 10px 2px 0;"><b>Path std:</b> {std_str}</td>
    </tr>
    <tr>
      <td style="padding:2px 10px 2px 0;"><b>LogLike:</b> {ll_str}</td>
      <td style="padding:2px 10px 2px 0;"><b>Multiplicity:</b> {multiplicity}</td>
      <td style="padding:2px 10px 2px 0;" colspan="2">
        <b>Nt values:</b> {nt_values}</td>
    </tr>
  </table>

  <div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:8px;">
    <div>
      <div style="font-size:11px;color:#888;margin-bottom:3px;text-align:center;">
        Viterbi spectrogram</div>
      {spec_img}
    </div>
    <div>
      <div style="font-size:11px;color:#888;margin-bottom:3px;text-align:center;">
        Fold &mdash; polynomial</div>
      {poly_img}
    </div>
    <div>
      <div style="font-size:11px;color:#888;margin-bottom:3px;text-align:center;">
        Fold &mdash; Kepler</div>
      {kep_img}
    </div>
  </div>

</div>""".format(
            cand_id=cand_id, card_border=card_border, card_bg=card_bg,
            rank=i+1, sel_badge=sel_badge, known_badge=known_badge,
            freq_str=freq_str, snr_str=snr_str, snr_poly=snr_poly,
            snr_kepler=snr_kepler, std_str=std_str, ll_str=ll_str,
            multiplicity=multiplicity, nt_values=nt_values,
            spec_img=img_tag(spec_path),
            poly_img=img_tag(poly_path),
            kep_img=img_tag(kep_path),
        )
        cards.append(card)

    # Index bar
    index_links = []
    for i, row in enumerate(candidates):
        cand_id   = row.get("cand_id", "")
        known     = row.get("known_match", "")
        snr       = row.get("snr_best", "")
        selection = row.get("selection", "")
        try:
            snr_f = "{:.1f}".format(float(snr))
        except ValueError:
            snr_f = snr
        col = ("#c8a800" if known
               else ("#81c995" if "binary" in selection else "#ce93d8"))
        index_links.append(
            '<a href="#{cid}" style="color:{col};font-size:11px;'
            'text-decoration:none;white-space:nowrap;">'
            '{rank}.{cid}({snr})</a>'.format(
                cid=cand_id, col=col, rank=i+1, snr=snr_f)
        )
    index_html = " &nbsp;".join(index_links)

    html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>47 Tuc Viterbi Candidate Gallery</title>
<style>
  body {{
    background:#111;color:#ddd;
    font-family:'Courier New',monospace;
    margin:0 auto;padding:20px;max-width:1600px;
  }}
  h1 {{color:#fff;border-bottom:1px solid #444;padding-bottom:8px;}}
  h2 {{color:#aaa;font-size:14px;font-weight:normal;}}
  .summary {{background:#1e1e1e;border:1px solid #333;border-radius:6px;
             padding:12px 16px;margin-bottom:16px;font-size:13px;}}
  .index   {{background:#161616;border:1px solid #333;border-radius:6px;
             padding:10px 14px;margin-bottom:24px;line-height:2.2;}}
  .legend span {{display:inline-block;padding:2px 8px;border-radius:3px;
                 font-size:12px;margin-right:8px;}}
</style>
</head>
<body>

<h1>47 Tuc Viterbi Blind Search &mdash; Candidate Gallery</h1>
<h2>47Tuc_blind_search_v1 &nbsp;|&nbsp; Sorted by S/N descending</h2>

<div class="summary">
  <b>Total candidates:</b> {n_total} &nbsp;|&nbsp;
  <b>Binary:</b> {n_binary} &nbsp;|&nbsp;
  <b>Isolated:</b> {n_isolated} &nbsp;|&nbsp;
  <b>Known pulsar matches:</b> {n_known}
</div>

<div class="legend" style="margin-bottom:14px;">
  Legend:
  <span style="background:#81c995;color:#000;">binary</span>
  <span style="background:#ce93d8;color:#000;">isolated</span>
  <span style="background:#4fc3f7;color:#000;">binary+isolated</span>
  <span style="background:#c8a800;color:#000;">known pulsar</span>
</div>

<div class="index">
<b style="color:#aaa;font-size:12px;">Jump to candidate:</b><br>
{index_html}
</div>

{cards}

</body>
</html>""".format(
        n_total=n_total, n_binary=n_binary,
        n_isolated=n_isolated, n_known=n_known,
        index_html=index_html,
        cards="".join(cards),
    )

    with open(args.out_html, "w") as fh:
        fh.write(html)

    size_mb = os.path.getsize(args.out_html) / 1e6
    print("\nGallery written: {}  ({:.1f} MB)".format(args.out_html, size_mb),
          flush=True)
    print("  Missing spectrograms: {}".format(n_missing_spec), flush=True)
    print("  Missing poly folds:   {}".format(n_missing_poly), flush=True)
    print("  Missing kepler folds: {}".format(n_missing_kep),  flush=True)
    print("\nDownload and open in any browser.", flush=True)


if __name__ == "__main__":
    main()