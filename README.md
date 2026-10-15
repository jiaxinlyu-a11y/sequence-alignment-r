# sequence-alignment-r

Pairwise sequence alignment in base R: **Needleman-Wunsch** (global) and **Smith-Waterman** (local). No packages needed. Each function fills the similarity matrix, traces back recursively, and returns **all** optimal alignments, not just one.

## Run

```
Rscript sequence-alignment.R
```

or `source("sequence-alignment.R")` in an R session. The script runs both algorithms on an example pair and prints the matrices, scores and alignments.

## Example

Sequences `s = AAGTGCCTCAAGATA` and `t = ACCGTCTCAGCAATA`, scoring: match +1, mismatch -1, gap -2.

| | Best score | Optimal alignments |
|---|---|---|
| Global (Needleman-Wunsch) | 0 | 16 (gaps in different places) |
| Local (Smith-Waterman) | 4 | 1: `CTCA` / `CTCA` |

One of the 16 global alignments:

```
AAGTGCCTCA-AGATA
-ACCGTCTCAGCAATA
```

## Global vs local

| | Global | Local |
|---|---|---|
| First row and column | 0, -2, -4, … | all 0 |
| Cell value | `max(diag, up, left)` | `max(0, diag, up, left)` |
| Best score | bottom-right cell | largest cell in the matrix |
| Traceback | from bottom-right to top-left | from each largest cell until a 0 |

## Functions

```r
needleman_wunsch(s, t, match = 1, mismatch = -1, gap = -2)
smith_waterman(s, t, match = 1, mismatch = -1, gap = -2)
```

Both take two strings and return a list with `matrix` (the similarity matrix), `score` (the best score) and `alignments` (a list; each item is a vector of two strings, the aligned `s` and `t`).

## Limits

Linear gap penalty only. No affine gaps and no substitution matrices such as BLOSUM, so it is meant for learning and small examples, not for large biological sequences.

## Files

- `sequence-alignment.R`: the code.
- `NOTES.md`: a step-by-step walkthrough of how the code works.
