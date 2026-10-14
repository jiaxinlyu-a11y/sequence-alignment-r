# Sequence alignment in base R: Needleman-Wunsch (global) and Smith-Waterman (local).
#
# Part 1: Needleman-Wunsch global alignment.
# Scoring: match = +1, mismatch = -1, gap = -2 (linear gap penalty).
# Builds the similarity matrix, then traces back ALL optimal alignments.

needleman_wunsch <- function(s, t, match = 1, mismatch = -1, gap = -2) {
  a <- strsplit(s, "")[[1]]
  b <- strsplit(t, "")[[1]]
  n <- length(a)
  m <- length(b)

  # Similarity matrix F: rows follow s, columns follow t, first row/column = gaps only
  F <- matrix(0, nrow = n + 1, ncol = m + 1)
  F[, 1] <- gap * (0:n)
  F[1, ] <- gap * (0:m)

  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      sub  <- if (a[i - 1] == b[j - 1]) match else mismatch
      diag <- F[i - 1, j - 1] + sub   # align a[i-1] with b[j-1]
      up   <- F[i - 1, j] + gap       # gap in t
      left <- F[i, j - 1] + gap       # gap in s
      F[i, j] <- max(diag, up, left)
    }
  }

  # Traceback: follow every move that produced F[i, j]; keep all optimal paths
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

  list(matrix = F, score = F[n + 1, m + 1], alignments = results)
}

# ---- Assignment 1 example ----
s <- "AAGTGCCTCAAGATA"
t <- "ACCGTCTCAGCAATA"

res <- needleman_wunsch(s, t)

cat("Global similarity matrix (rows: s, columns: t)\n")
print(res$matrix)
cat("\nOptimal global score:", res$score, "\n")
cat("Number of optimal alignments:", length(res$alignments), "\n\n")
for (k in seq_along(res$alignments)) {
  cat("Alignment", k, "\n")
  cat(res$alignments[[k]][1], "\n")
  cat(res$alignments[[k]][2], "\n\n")
}


# =============================================================================
# Smith-Waterman local alignment in base R.
# Same scoring as above. Differences from Needleman-Wunsch:
#   1. first row and column are 0 (leading gaps cost nothing)
#   2. every cell is at least 0 (a bad stretch is dropped and we start over)
#   3. the answer is the largest cell, and the traceback stops at a 0
# =============================================================================

smith_waterman <- function(s, t, match = 1, mismatch = -1, gap = -2) {
  a <- strsplit(s, "")[[1]]
  b <- strsplit(t, "")[[1]]
  n <- length(a)
  m <- length(b)

  # First row and column stay 0
  F <- matrix(0, nrow = n + 1, ncol = m + 1)

  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      sub  <- if (a[i - 1] == b[j - 1]) match else mismatch
      diag <- F[i - 1, j - 1] + sub
      up   <- F[i - 1, j] + gap
      left <- F[i, j - 1] + gap
      F[i, j] <- max(0, diag, up, left)   # the extra 0 is what makes it local
    }
  }

  best   <- max(F)                          # best local score, anywhere in F
  starts <- which(F == best, arr.ind = TRUE) # every cell that reaches it

  # Traceback: same three moves as before, but stop when the score is 0
  results <- list()
  trace <- function(i, j, top, bottom) {
    if (F[i, j] == 0) {
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
  for (r in seq_len(nrow(starts))) {
    trace(starts[r, 1], starts[r, 2], "", "")
  }

  list(matrix = F, score = best, alignments = results)
}

# ---- Same sequences, local alignment ----
res_local <- smith_waterman(s, t)

cat("\n\nLocal similarity matrix (rows: s, columns: t)\n")
print(res_local$matrix)
cat("\nBest local score:", res_local$score, "\n")
cat("Number of optimal local alignments:", length(res_local$alignments), "\n\n")
for (k in seq_along(res_local$alignments)) {
  cat("Local alignment", k, "\n")
  cat(res_local$alignments[[k]][1], "\n")
  cat(res_local$alignments[[k]][2], "\n\n")
}
