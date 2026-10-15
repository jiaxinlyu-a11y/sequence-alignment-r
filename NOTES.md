# Sequence alignment notes: Needleman-Wunsch and Smith-Waterman (in R)

Code: `sequence-alignment.R` (one file with a global and a local alignment function)
Example sequences: `s = AAGTGCCTCAAGATA`, `t = ACCGTCTCAGCAATA`
Scoring: match +1, mismatch -1, gap -2

---

## 1. What the two algorithms do

| | Global: Needleman-Wunsch | Local: Smith-Waterman |
|---|---|---|
| Goal | Align both sequences end to end | Find the most similar stretch only |
| First row and column | 0, -2, -4, -6 … (leading gaps cost points) | All 0 |
| Lowest value in a cell | Can be negative | 0 |
| Best score is | The bottom-right cell | The largest cell in the matrix |
| Traceback starts at | Bottom-right cell | Each cell holding the maximum (there may be several) |
| Traceback stops at | Top-left cell | A cell with value 0 |

Results on the example:

- Global: best score **0**, with **16** equally good alignments (the gaps sit in different places).
- Local: best score **4**, **1** alignment: `CTCA` against `CTCA`.
- The global matrix matches my assignment result cell by cell, and both alignments listed in the assignment are among the 16.

---

## 2. Code walkthrough (Needleman-Wunsch)

### 2.1 The start of the function

```r
needleman_wunsch <- function(s, t, match = 1, mismatch = -1, gap = -2) {
  a <- strsplit(s, "")[[1]]
  b <- strsplit(t, "")[[1]]
  n <- length(a)
  m <- length(b)
```

- The last three arguments are the scoring rules. The defaults are the assignment's +1, -1, -2; to change the scoring, change them here.
- `strsplit(s, "")` splits on the empty string, i.e. **one character at a time**.
- `[[1]]`: `strsplit` returns a **list** (it can split several strings at once). We split one string, so we take the first item of the list, which is an ordinary character vector. `[[1]]` has nothing to do with "row 1, column 1" of a matrix.

```r
strsplit("AAG", "")
# [[1]]
# [1] "A" "A" "G"
```

- A **string** is not the same as a **character vector**:

```r
s <- "AAGTGCCTCAAGATA"
length(s)    # 1: the vector holds one item (the whole string)
nchar(s)     # 15: the string has 15 characters
a <- strsplit(s, "")[[1]]
length(a)    # 15: the vector holds 15 items, one character each
a[3]         # "G"
```

- Splitting into a vector lets us write `a[i - 1] == b[j - 1]` and compare one character at a time.
- `n` and `m` are the sequence lengths, 15 here.

### 2.2 Building and initialising the matrix

```r
F <- matrix(0, nrow = n + 1, ncol = m + 1)
F[, 1] <- gap * (0:n)
F[1, ] <- gap * (0:m)
```

- The matrix is **(n+1) × (m+1)**. The extra first row and column stand for the empty sequence.
- The first row and column are **not all zeros**; only the top-left cell is 0. A sequence against an empty sequence can only be gaps, so each extra character costs one gap: 0, -2, -4, -6, …
- `F[, 1]` is column 1 (leaving the part before the comma empty means all rows). `0:n` makes the integers 0 to n. `gap * (0:n)` is 0, -2, -4, …, -30.
- Note: **a first row and column of all zeros belongs to local alignment (Smith-Waterman)**. Don't mix the two up.

### 2.3 Filling the matrix (the core)

```r
for (i in 2:(n + 1)) {
  for (j in 2:(m + 1)) {
    sub  <- if (a[i - 1] == b[j - 1]) match else mismatch
    diag <- F[i - 1, j - 1] + sub
    up   <- F[i - 1, j] + gap
    left <- F[i, j - 1] + gap
    F[i, j] <- max(diag, up, left)
  }
}
```

- We start at row 2, column 2, because the first row and column are already filled. Each cell depends only on the cells to its **upper left, directly above, and directly left**, all of which are already computed.
- Why `i - 1` and `j - 1` for the characters: the matrix has an extra row and column, so row `i` of the matrix corresponds to the `(i-1)`th character of `a`.
- Meaning of `F[i, j]`: the best score for globally aligning the first i-1 characters of s with the first j-1 characters of t. The bottom-right cell is the best score for the whole pair.

**The three moves** (rule of thumb: **only the diagonal looks at whether the characters match; the other two directions just add the gap**):

| Move | Meaning | Score added |
|---|---|---|
| From the diagonal (diag) | The two characters are paired | +1 if equal, -1 if not |
| From above (up) | A character of s faces a gap | gap (-2) |
| From the left (left) | A character of t faces a gap | gap (-2) |

**Why up and left only add the gap:** picture the alignment as two rows, scored column by column. A column with two characters is scored by whether they match. A column with one character facing `-` has nothing to compare, so it is simply charged the gap. Every step in the matrix adds one column to the alignment.

**Picturing it in the matrix:** **rows are s, columns are t.** Imagine two fingers, one on s and one on t:

| Move in the matrix | New column in the alignment (s on top, t below) |
|---|---|
| Diagonal | character of s / character of t |
| Down | character of s / `-` (s moved, t did not) |
| Right | `-` / character of t (t moved, s did not) |

The sequence whose finger moved gets a character in that column; the one that did not move gets `-`.

**A small example:** `s = "AG"`, `t = "A"`

```
          empty  A
  empty     0   -2
  A        -2    ?
  G        -4    ?
```

- `F[2,2]`: diag = 0 + 1 = 1; up = -2 - 2 = -4; left = -2 - 2 = -4. The maximum is **1**.
- `F[3,2]`: diag = `F[2,1]` + (-1) = -3; up = `F[2,2]` + gap = 1 - 2 = -1; left = `F[3,1]` + gap = -6. The maximum is **-1** (from up).
- The matching alignment is `AG` over `A-`, scoring +1 - 2 = -1, which agrees with the matrix.

### 2.4 Traceback: recovering the alignment

Once the matrix is filled, walk backwards from the bottom-right cell to the top-left. At each cell ask: **how did this number come about?** Whichever of diagonal, above or left produced it, step back that way and add one column to the alignment.

```r
results <- list()
trace <- function(i, j, top, bottom) {
  if (i == 1 && j == 1) {
    results[[length(results) + 1]] <<- c(top, bottom)
    return(invisible())
  }
  if (i > 1 && j > 1) {
    sub <- if (a[i - 1] == b[j - 1]) match else mismatch
    if (F[i, j] == F[i - 1, j - 1] + sub)
      trace(i - 1, j - 1, paste0(a[i - 1], top), paste0(b[j - 1], bottom))
  }
  if (i > 1 && F[i, j] == F[i - 1, j] + gap)
    trace(i - 1, j, paste0(a[i - 1], top), paste0("-", bottom))
  if (j > 1 && F[i, j] == F[i, j - 1] + gap)
    trace(i, j - 1, paste0("-", top), paste0(b[j - 1], bottom))
}
trace(n + 1, m + 1, "", "")
```

Point by point:

- `trace` is a **recursive function**: after each step back it calls itself again from the new position, until it reaches the top-left. `i, j` is the current cell; `top` and `bottom` are the alignment built so far (s on the upper row, t on the lower row).
- The last line starts at the bottom-right cell `(n+1, m+1)` with both rows empty strings `""`.
- `if (i == 1 && j == 1)`: we reached the top-left, so the path is complete. Store the alignment in `results` and stop.
- `&&` means "and" (logical AND): true only if both sides are true. Use it inside `if`. A single `&` also means "and" but works element by element on whole vectors. "Or" is `||` and `|` in the same way.
- `results[[length(results) + 1]] <<- c(top, bottom)`:
  - `length(results) + 1` is the number of the next free slot. This **stores** something in that slot; it does not read anything.
  - `c(top, bottom)` combines the two alignment rows into a vector of 2 items. **What is stored is the two aligned sequences, not a score.**
  - `<<-` writes to `results` outside the function. A plain `<-` would only change a temporary copy inside the function.
- `return(invisible())`: ends the function and prints nothing. The point is to stop once the top-left is reached.
- `paste0(a[i - 1], top)`: joins a character and a whole piece of text end to end into a new string. Because we **walk backwards**, each new character goes at the **front**, so the finished alignment reads in the right order. Writing it the other way round (`paste0(top, a[i - 1])`) would give a reversed alignment.
- The three `if` statements are **independent**, not `else if`. If more than one direction holds for a cell, every one of them is followed and each continues on its own. This is why there are **16 optimal alignments**.

**Key point: why are `top` and `bottom` strings?**
No type is declared anywhere. It is decided by the initial values passed in the call `trace(n + 1, m + 1, "", "")`: two `""`, i.e. empty strings, so they stay strings all the way. R needs no type declarations; a variable is whatever you give it. Calling `trace(n + 1, m + 1, c(), c())` would make them empty vectors instead, and the later code would have to change accordingly.

**What up and left mean:**

- Coming from above: s used one more character, and t gets `-` in that column (scored with the gap, not by matching a character of t).
- Coming from the left: t used one more character, and s gets `-` in that column.

### 2.5 Returning the result

```r
list(matrix = F, score = F[n + 1, m + 1], alignments = results)
}
```

- In R the last expression of a function is its return value, so `return` is not needed.
- A function can hand back only one object, so the three results are packed into a list and read back by name with `$`:

```r
res <- needleman_wunsch(s, t)
res$matrix                  # the matrix
res$score                   # best score
res$alignments              # all alignments
res$alignments[[1]]         # alignment 1 (a vector of 2 items)
res$alignments[[1]][1]      # alignment 1, the s row
res$alignments[[1]][2]      # alignment 1, the t row
```

- `length(res$alignments)` is **how many items the list has** (16 here). It is not a table of "16 rows and 2 columns": each item of the list is a vector of length 2. A list has items, not rows and columns.
- In global alignment the bottom-right cell is always the best score for the whole pair.

### 2.6 The printing part

```r
for (k in seq_along(res$alignments)) {
  cat("Alignment", k, "\n")
  cat(res$alignments[[k]][1], "\n")
  cat(res$alignments[[k]][2], "\n\n")
}
```

- `cat` prints text and `\n` is a line break. `seq_along(x)` gives `1, 2, ..., length(x)`, and `k` takes each value in turn.

---

## 3. Smith-Waterman (local): three differences from global

```r
F <- matrix(0, nrow = n + 1, ncol = m + 1)   # 1. first row and column simply stay 0
F[i, j] <- max(0, diag, up, left)             # 2. an extra 0: drop a bad stretch and start over
best   <- max(F)                              # 3. the answer is the largest cell in the matrix
starts <- which(F == best, arr.ind = TRUE)
```

- The local matrix has **no negative values**.
- The traceback **stops at a 0** (global stops only at the top-left).

### `which(F == best, arr.ind = TRUE)` and the `for` loop

```r
best   <- max(F)
starts <- which(F == best, arr.ind = TRUE)
for (r in seq_len(nrow(starts))) {
  trace(starts[r, 1], starts[r, 2], "", "")
}
```

- `F == best` compares every cell with the maximum and gives a TRUE/FALSE matrix of the same size.
- `which(..., arr.ind = TRUE)` finds the TRUE positions and writes them as row and column numbers.

```r
M <- matrix(c(0, 1, 4, 2, 4, 0), nrow = 2)
which(M == 4)                    # [1] 3 5    without arr.ind: just position numbers, awkward to use
which(M == 4, arr.ind = TRUE)    # with it: a matrix with columns row and col
#      row col
# [1,]   1   2
# [2,]   1   3
```

- **`starts` is a matrix, not a list.** That is why we read it with `starts[r, 1]` (row number of the r-th start) and `starts[r, 2]` (column number); lists use `[[ ]]`.
- What makes it a matrix is `arr.ind = TRUE`. In R, **a variable's type is decided by the result assigned to it**.
- The maximum can occur several times, and each occurrence ends one optimal local alignment, so the `for` loop runs a traceback from every such cell.
- In our example the maximum, 4, occurs once, at row 11, column 10. The loop runs once, equivalent to `trace(11, 10, "", "")`, which steps diagonally back to a 0 and gives `CTCA`.

---

## 4. R reminders

- `strsplit(x, "")[[1]]` turns a string into a character vector. `length()` counts items in a vector; `nchar()` counts characters in a string.
- `&&` and `||` are for single TRUE/FALSE values (use them in `if`); `&` and `|` work element by element on vectors.
- `paste0(a, b)` joins end to end with no space.
- `c(x, y)` strings things together into a vector.
- `<<-` changes a variable outside the function.
- `return(invisible())` ends a function without printing.
- Vector, matrix, list: a vector is a row of items of the same type; a matrix is a regular table read with `[i, j]`; a list is a series of slots that can hold different things, read with `[[k]]`.
- A variable's type comes from the value assigned (or passed in); no declaration is needed.

### Getting R output into a PDF

- Put the code in a `.qmd` file with `format: pdf` in the header (use `format: typst` if LaTeX is not installed) and click Render. You get a PDF with the code and its output.
- A plain `.R` file also works: `rmarkdown::render("file.R", output_format = "pdf_document")`, or Ctrl+Shift+K in RStudio. For text to appear in the body of the PDF, write comments as `#'`.

---

## 5. Talking points for an interview

- I implemented global and local alignment from scratch in R: fill the similarity matrix, then trace back recursively and list **all** optimal alignments rather than just one.
- Global and local differ in only three places: the initialisation, taking `max(0, ...)` in each cell, and taking the maximum cell as the answer and stopping at 0.
- I checked it against my coursework sequences: the global matrix matches the assignment cell by cell, the best score is 0 with 16 equally good alignments; the local score is 4 and the alignment is `CTCA`.
