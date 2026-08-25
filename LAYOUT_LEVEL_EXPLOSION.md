# Why the fast layout produces far too many lower-level IBFs

Analysis of `generate_layout` at commit `adb502a` (source unchanged since `4fa995e`).

## TL;DR

The number of `#LOWER_LEVEL_IBF` lines explodes because level-1 IBFs fail to come out
as pure *split* layouts. Four defects, in decreasing order of impact:

1. **The overflow sentinel becomes a layout.** When `refine_and_bin` never finds a feasible
   `t_max`, `binning_core` returns `{0, 0, 0}` as `(split_start, split_bins, merge_start)`.
   With `merge_start == 0`, `generate_hibf` turns *every* non-empty technical bin into its own
   lower-level IBF — including single-sequence bins. One such IBF costs ~`global_bins` spurious
   children. This is almost certainly what produces the 5901.
2. **`refine_and_bin` returns the last `t_max` it tried, not the best one.** The `p`-th result
   is returned without ever being scored, and `b_res` is only consulted when the final result
   overflowed.
3. **The scoring metric actively drives `t_max` onto the cliff** where merge bins first appear,
   because `merge_average()` returns `0` when there are no merge bins and the comparison reads
   that as "`t_max` too small".
4. **`global_bins` is fixed for the whole tree**, so the root can pack more sequences into a
   merge bin than the child IBF has technical bins, forcing merging again one level down.

Defects 2 + 3 together can produce an **exact fixed point**: an IBF whose only merge bin
contains precisely the sequence set it was given. `generate_hibf` then recurses on that same
set forever. **On the 1024 dataset the current code does not terminate.**

## Reproduction

```
cd executables && make -B
./generate_layout /srv/data/smehringer/data_simulated/hibf_paper_lsh/$BINS \
    8 19 19 6767 0.05 15:4,20:5,15:6,10:12 0.001 5 2 $BINS.output.layout 32
```

Counting with `grep -c LOWER_LEVEL_IBF`:

| `$BINS` | observed | expected (chopper) |
|--------:|---------:|-------------------:|
| 1024    | 64       | 64  |
| 2048    | 64       | 64  |
| 4096    | 64       | 64  |
| 8192    | 128      | 128 |
| 16384   | 5901     | 128 |
| 32768   | 1178     | 192 |
| 65536   | 3014     | 256 |
| 131072  | 4426     | 384 |
| 262144  | 8105     | 512 |
| 524288  | 15277    | 768 |

Measured here with an instrumented build (`-O2`, per-IBF logging of `n_seqs`, `t_max`,
split/merge bin counts and child counts):

* **1024** — does not terminate. Still emitting 12 identical IBFs per level at level 21
  when killed. Reproduced identically at `9728f90`, so this predates the merge/seeding reorder.
* **16384** — terminates with **604** lower-level IBFs (128 at level 1 + 476 at level 2).

604 vs. your 5901 is run-to-run variance: `binning_core` iterates `std::unordered_map`s
(`labMaps`), so bin composition depends on hash/allocator order, which differs between
`-O2` and your `-O3 -flto` build. The mechanism is the same; only how many IBFs land in the
bad regime changes. See "Where 5901 comes from" below.

## Background: what this dataset looks like

`hibf_paper_lsh/$BINS` is always the same ~65 GB of content, cut into `$BINS/64` genome groups
of 64 pieces each:

| `$BINS` | groups | avg file size |
|--------:|-------:|--------------:|
| 1024    | 16     | 63.3 MB |
| 4096    | 64     | 15.8 MB |
| 16384   | 256    |  3.96 MB |
| 65536   | 1024   |  0.99 MB |

Two consequences:

* The 64 files within a group are **near-duplicates** — a whole group costs about as much
  as a single file. Confirmed by the sketch sizes: for `1024`, union = 996 118 vs.
  sum = 63 983 552 (ratio 1:64).
* LSH finds exactly one cluster per group, at *every* level:
  `[lsh] levels: lvl0=16 lvl1=16 lvl2=16 lvl3=16` for 1024,
  `lvl0=256 ... lvl3=256` for 16384. The hierarchy is degenerate on this data — expected,
  since the 64 cuts genuinely are near-identical, but it means `binning_core` effectively
  works with one flat clustering and the "climbing" mechanism has nowhere to descend.

The target shape is two levels: a root with `global_bins` merge bins, each child holding
~`global_bins` sequences that all become split bins. For 16384 that is exactly
128 root bins × 128 sequences → 128 lower-level IBFs, no level 2.

---

## Cause 1 — the overflow sentinel is returned as a real layout

`code/templates/fast_construct_binning_core.tpp:140` and
`code/templates/fast_construct_binning_core.tpp:182` both return

```cpp
return {res, {0,0,0}, track_fill, true};
```

on overflow. `code/templates/fast_construct_generate_hibf.tpp:132` only falls back to the last
good result when one exists:

```cpp
auto final_res = (final_is_overflow && valid) ? b_res : res;
```

If `valid` was never set — i.e. *every* `t_max` the search could reach overflowed — the
sentinel is handed back as the layout. `code/templates/fast_construct_generate_hibf.tpp:213`
then walks from `merge_start == 0`:

```cpp
for (size_t b = merge_start; b < ibf.size(); b++) {
    if (ibf[b].empty()) continue;
    tasks.push_back({ibf_index, b});
}
```

so every non-empty technical bin — including bins holding a single sequence — becomes a
lower-level IBF.

**Why the search can fail entirely.** The ceiling is set at
`code/templates/fast_construct_generate_hibf.tpp:48`:

```cpp
curr_upper = curr_t_max * 2;     // == 2 * (union + sum) / (2 * bins * s)
```

That ceiling lies *below the largest single user bin* whenever an IBF has fewer sequences
than `bins`. Every reachable `t_max` then force-splits every sequence at
`code/templates/fast_construct_binning_core.tpp:142` (`est_sketch_size >= t_max`), the split
cascade exhausts the bins, and overflow is unavoidable:

```
[L1#0] START n_seqs=67 sub_bins=128 lower=15638 upper=272175 init_t_max=1124269
[L1#0]   OVERFLOW t_max=1124269
[L1#0]   OVERFLOW t_max=1686403
...                                   (22 retries, converging on the ceiling)
[L1#0]   OVERFLOW t_max=2248537
[L1#0] END overflow=1 merge_start=0 -> children=128
```

Per-sequence size here is ~4.06e6 while the reachable ceiling is 2.25e6, so all 67 sequences
split into ≥2 bins → 134 > 128 → overflow at every step. One IBF, 128 spurious children.

This path is also a **correctness bug**, not only a size bug: the sentinel `res` is the
*partial* binning built before the loop bailed out, so sequences that were never placed are
silently dropped from the layout.

## Cause 2 — `refine_and_bin` returns the last `t_max`, not the best

`code/templates/fast_construct_generate_hibf.tpp:89` breaks out **before** the result is scored:

```cpp
if(it == p) break;          // res at the (p+1)-th t_max is never evaluated
```

and `b_res` / `best_t_max`, set at
`code/templates/fast_construct_generate_hibf.tpp:94-96`, are only used on overflow
(`code/templates/fast_construct_generate_hibf.tpp:132-133`). So the returned layout is the one
binning nobody looked at.

All 12 affected level-1 IBFs in the 1024 run look identical:

```
[L1#15]   [it 4] t_max=62122583 used=64 split_bins=64 merge_bins=0   <- ideal, sits in b_res
[L1#15] END      t_max=62615619 nonempty=1 merge_start=63 -> children=1
```

Iteration 4 produced the perfect answer (64 split bins, zero children). The 6th, unscored
evaluation at a `t_max` 0.8 % higher collapsed all 64 sequences into one bin — and that is
what was returned.

## Cause 3 — the metric steers `t_max` straight onto the cliff

`code/templates/fast_construct_generate_hibf.tpp:119`:

```cpp
if (split_bin_amt == 0 || split_avg < merge_avg) curr_upper = curr_t_max;  // lower t_max
else                                             curr_lower = curr_t_max;  // raise t_max
```

When an IBF has **no merge bins** — the outcome we want — `merge_average` returns `0`
(`code/fast_construct_bin.cpp:24`), so `split_avg < merge_avg` is always false and `t_max` is
raised. The search is therefore guaranteed to climb until merge bins appear: it converges on
the smallest `t_max` that creates children.

For near-identical sequences that boundary is a **cliff, not a gradient**:

* just below it, every sequence is force-split → all split bins, zero children;
* just above it, `t_max > cluster union size`, and the following cascade fires inside
  `binning_core`:
  1. `code/templates/fast_construct_binning_core.tpp:255-259` — the whole cluster is classified
     "small", so all its sequences go into `reserved_small_seqs`;
  2. `code/templates/fast_construct_binning_core.tpp:265` —
     `allowed_merge = union(supersketch) / (s * t_max)` truncates to **0**, so
     `seeding_max_bin == bins` (`:271`) and the entire late-merge round-robin at `:275` is
     skipped;
  3. `code/templates/fast_construct_binning_core.tpp:315-318` — every bottom-level cluster is
     `all_binned` (everything is "reserved"), so there are no candidates; seeding (`:361-363`)
     and climbing (`:403`) place nothing;
  4. the fallback at `code/templates/fast_construct_binning_core.tpp:438` dumps everything into
     one bin. Near-duplicates add ~0 new k-mers, so the `t_max` guard in
     `try_insert_sequence` (`code/templates/fast_construct_binning_core.tpp:58`) never fires,
     and after the first insertion `bin_for_cluster` funnels the rest into the same bin.

**That bin then contains exactly the sequence set the IBF was given**, and `generate_hibf`
recurses on it unconditionally. It is an exact fixed point:

```
[L1#15] START n_seqs=64 sub_bins=64 lower=62122 upper=3976832 init_t_max=31554328
[L1#15] END   ... nonempty=1 merge_start=63 -> children=1
[L2#0]  START n_seqs=64 sub_bins=64 lower=62122 upper=3976832 init_t_max=31554328   <- identical
[L2#0]  END   ... nonempty=1 merge_start=63 -> children=1
[L3#0]  START n_seqs=64 sub_bins=64 lower=62122 upper=3976832 init_t_max=31554328   <- identical
```

Level counts on 1024:

```
[LEVEL 1] ibfs=27 -> tasks=12 total_seqs_in_tasks=768
[LEVEL 2] ibfs=12 -> tasks=12 total_seqs_in_tasks=768
...
[LEVEL 21] ibfs=12 -> tasks=12 total_seqs_in_tasks=768     <- still going when killed
```

Whether this fixed point is hit or narrowly escaped depends on the ratio `bins / n_seqs`:

* 1024, child with 64 seqs in 64 bins: ceiling `2 * init = 63 108 656` vs. per-sequence size
  `62 138 000` → ceiling is *above* → no split → collapse into one merge bin → infinite.
* 16384, child with 64 seqs in 128 bins: ceiling `1 911 914` vs. per-sequence size
  `3 765 000` → ceiling is *below* → each sequence splits in two → 128 split bins, zero
  children → terminates.

Termination is currently a coin flip.

## Cause 4 — one global bin count for the whole tree

`code/templates/fast_construct_generate_hibf.tpp:168-176` computes `global_bins` once from the
total number of user bins and reuses it for every IBF at every level
(`code/templates/fast_construct_generate_hibf.tpp:193` and `:244`).

Nothing then constrains a merge bin to hold at most `global_bins` sequences. On 16384 the root
packs 130–131 sequences into bins whose child IBF has only 128 technical bins:

```
n_seqs per level-1 IBF:  67, 84, 89, 90, 95, 111, 113, 116, then 130 (x101), 131 (x19)
```

Each of the 120 over-full children *must* merge again, contributing 1–5 level-2 IBFs
(348 of them in my run). The eight under-full ones are exactly the ones at risk of Cause 1.

---

## Where 5901 comes from

Breakdown of my 16384 run (604 lower-level IBFs):

| source | count |
|---|---:|
| level 1, legitimate | 128 |
| level 2, from the one IBF that returned the overflow sentinel (Cause 1) | 128 |
| level 2, from the 120 over-packed children (Cause 4) | 348 |

Scaling that to your number: `128 * k + ~140 ≈ 5901` gives `k ≈ 45`, i.e. roughly 45 of your
128 level-1 IBFs returned the overflow sentinel instead of 1 of mine. That is consistent with
`-O3 -flto` producing a different hash-map iteration order and hence a different distribution
of sequences over the root's bins. **Cause 1 dominates your numbers.**

## Verification

Patching only two things in a scratchpad copy of the source:

* **A** — score every feasible result in the refinement loop (including the `p`-th) and return
  the best, instead of returning the last;
* **B** — raise `curr_upper` at `code/templates/fast_construct_generate_hibf.tpp:48` above the
  largest single user bin, so a feasible `t_max` always exists:
  `curr_upper = max(curr_t_max * 2, max_i(|sketch_i| / s) * 17/16 + 1)`

gives:

| dataset | current | patched |
|---|---|---|
| 1024 | never terminates | terminates |
| 16384 | 604 lower-level IBFs, 1 overflow sentinel | 367, **0** overflow sentinels |

367 is still far from the ideal 128 — the throwaway score in patch A merely minimizes the
child count, so it picks a huge `t_max` and uses 9 of 64 root bins, and Cause 4 is untouched.
The experiment confirms the mechanisms; it is not a usable fix.

## Suggested fixes

1. **Make overflow a failure, not a layout.** `binning_core` should not return `{0,0,0}` as a
   usable range tuple. If `refine_and_bin` cannot find any feasible `t_max`, that is a bug to
   escalate (or retry with a raised ceiling), never something to emit. Right now it both
   inflates the level count and silently drops sequences.
2. **Fix the `t_max` search bounds** (`code/templates/fast_construct_generate_hibf.tpp:48`).
   The ceiling must be at least `max_i(|sketch_i| / s)`, so that "no user bin is force-split"
   is always reachable. The current `(union + sum) / (2 * bins * s)` scales with the IBF's
   content but `bins` does not, which is why small children can never find a feasible `t_max`.
3. **Give the refinement loop a real objective.** Score every candidate — bin size × number of
   children is the quantity the HIBF actually pays for — and return the best. Treat
   `merge_bin_amt == 0` as "done" rather than as "`t_max` too small"; the current
   `merge_avg == 0` comparison at
   `code/templates/fast_construct_generate_hibf.tpp:119` guarantees the search walks onto the
   collapse point.
4. **Add a progress guard in `generate_hibf`.** A child must never be handed the exact
   sequence set its parent bin held. If `binning` returns a single merge bin equal to its
   input, either force a split or stop recursing. This makes termination structural rather
   than accidental.
5. **Cap merge-bin occupancy** at the number of technical bins the child IBF will have, or
   size sub-IBFs from their own content instead of reusing `global_bins`
   (`code/templates/fast_construct_generate_hibf.tpp:176`) at every level.

## Instrumentation used

The measurements above come from a copy of `code/` with logging added to
`generate_hibf` — per `refine_and_bin` call: tag, `n_seqs`, `sub_bins`, `lower`, `upper`,
initial and chosen `t_max`, every overflow retry, per-iteration split/merge bin counts and
averages, and a final line with `nonempty`, `split_start`, `merge_start` and the resulting
child count; plus a per-level line with the number of IBFs and tasks. `generate_layout.cpp`
also gained a sketch cache (an optional 13th argument) so that runs after the first skip the
FASTA hashing entirely, which is what makes iterating on this practical.
