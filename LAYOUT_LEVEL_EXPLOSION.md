# Why the fast layout produced far too many lower-level IBFs

Diagnosis of `generate_layout` at commit `adb502a` (source unchanged since `4fa995e`), and the
fix, which is commit `a0af981`.

## Result

`grep -c LOWER_LEVEL_IBF` on the layout of each dataset:

| `$BINS` | before | after | chopper | levels after |
|--------:|-------:|------:|--------:|-------------:|
| 1024    | *never terminates* | 16   | 64  | 2 |
| 2048    | 64     | 32    | 64  | 2 |
| 4096    | 64     | 64    | 64  | 2 |
| 8192    | 128    | 128   | 128 | 2 |
| 16384   | 5901   | 128   | 128 | 2 |
| 32768   | 1178   | 192   | 192 | 2 |
| 65536   | 3014   | 256   | 256 | 2 |
| 131072  | 4426   | 384   | 384 | 2 |
| 262144  | 8105   | 514   | 512 | 3 |
| 524288  | 15277  | 785   | 768 | 3 |

Every layout was checked three ways: chopper's own `display_layout general` parses them (for
16384 it reports 128 top-level bins, all `merged`, `ub_count=128`, sizes 7.81e6 within 0.5% of
each other); a structural checker confirms that every user bin appears exactly once, that
header lines and content agree and are in level order, and that technical bin spans do not
overlap; and `function_tests` produces byte-identical output to a build from before the fix.

## The four defects

The number of `#LOWER_LEVEL_IBF` lines exploded because level-1 IBFs failed to come out as pure
*split* layouts. Four defects, in decreasing order of impact:

1. **The overflow sentinel becomes a layout.** When `refine_and_bin` never finds a feasible
   `t_max`, `binning_core` returns `{0, 0, 0}` as `(split_start, split_bins, merge_start)`.
   With `merge_start == 0`, `generate_hibf` turns *every* non-empty technical bin into its own
   lower-level IBF — including single-sequence bins. One such IBF costs ~`global_bins` spurious
   children. This is what produced the 5901.
2. **`refine_and_bin` returns the last `t_max` it tried, not the best one.** The `p`-th result
   is returned without ever being scored, and `b_res` is only consulted when the final result
   overflowed.
3. **The scoring metric actively drives `t_max` onto the cliff** where merge bins first appear,
   because `merge_average()` returns `0` when there are no merge bins and the comparison reads
   that as "`t_max` too small".
4. **`global_bins` is fixed for the whole tree**, so the root can pack more sequences into a
   merge bin than the child IBF has technical bins, forcing merging again one level down.

Defects 2 + 3 together produced an **exact fixed point**: an IBF whose only merge bin contained
precisely the sequence set it was given. `generate_hibf` then recursed on that same set forever.
**On the 1024 dataset the old code did not terminate.**

## Reproduction

```
cd executables && make -B
./generate_layout /srv/data/smehringer/data_simulated/hibf_paper_lsh/$BINS \
    8 19 19 6767 0.05 15:4,20:5,15:6,10:12 0.001 5 2 $BINS.output.layout 32
```

The binary now takes an optional 13th argument, a sketch cache file — see *Testing* at the end.

Counting with `grep -c LOWER_LEVEL_IBF`, **before the fix**:

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

604 against the 5901 reported from the original runs is run-to-run variance: `binning_core` iterates `std::unordered_map`s
(`labMaps`), so bin composition depends on hash/allocator order, which differs between
the `-O2` instrumented build and the `-O3 -flto` one. The mechanism is the same; only how many IBFs land in the
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

## Defect 1 — the overflow sentinel is returned as a real layout

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

## Defect 2 — `refine_and_bin` returns the last `t_max`, not the best

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

## Defect 3 — the metric steers `t_max` straight onto the cliff

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

## Defect 4 — one global bin count for the whole tree

`code/templates/fast_construct_generate_hibf.tpp:168-176` computes `global_bins` once from the
total number of user bins and reuses it for every IBF at every level
(`code/templates/fast_construct_generate_hibf.tpp:193` and `:244`).

Nothing then constrains a merge bin to hold at most `global_bins` sequences. On 16384 the root
packs 130–131 sequences into bins whose child IBF has only 128 technical bins:

```
n_seqs per level-1 IBF:  67, 84, 89, 90, 95, 111, 113, 116, then 130 (x101), 131 (x19)
```

Each of the 120 over-full children *must* merge again, contributing 1–5 level-2 IBFs
(348 of them in the instrumented run). The eight under-full ones are exactly the ones at risk
of defect 1.

---

## Where 5901 comes from

Breakdown of the instrumented 16384 run (604 lower-level IBFs):

| source | count |
|---|---:|
| level 1, legitimate | 128 |
| level 2, from the one IBF that returned the overflow sentinel (defect 1) | 128 |
| level 2, from the 120 over-packed children (defect 4) | 348 |

Scaling that up: `128 * k + ~140 ≈ 5901` gives `k ≈ 45`, i.e. roughly 45 of the 128 level-1
IBFs in the original run returned the overflow sentinel, against 1 here. That is consistent with
`-O3 -flto` producing a different hash-map iteration order and hence a different distribution
of sequences over the root's bins. **Defect 1 dominated the numbers.**

## The fix

All of this is commit `a0af981`. Line numbers below refer to that commit.

**Defect 2 and 3 — the refinement loop was rewritten** in
`code/templates/fast_construct_generate_hibf.tpp`:

* `layout_cost` (`:40`) scores a candidate layout in k-mers:
  `technical_bins * fullest_bin + sum of merge bin content`. The first term is what the IBF
  allocates — all of its bins at the size of the fullest — and the second is a lower bound on
  what the children below it will cost.
* `consider` (`:134`) scores *every* feasible candidate, including the `p`-th, and keeps the
  best. The result that ships is no longer the one nobody looked at.
* The direction rule (`:151`) is now simply "overflow means `t_max` is too small, raise it;
  otherwise lower it", because a smaller `t_max` makes every technical bin cheaper. The
  `merge_avg == 0` comparison that guaranteed a walk onto the collapse point is gone, and with
  it the dependency on `merge_average` / `splitting_average`.

**Defect 1 — overflow can no longer become a layout.** The search ceiling (`:118`, `:121`) is
now

```cpp
const size_t feasible_t_max = max_single + max_single/16 + 1;   // max_single = max_i(|sketch_i|/s)
size_t curr_upper = std::max(curr_t_max * 2, feasible_t_max);
```

Above `max_single` no user bin is split at all, so `binning_core` cannot run out of technical
bins: a non-overflowing `t_max` is always inside the interval. If the loop still finds nothing
usable, `refine_and_bin` runs once at `feasible_t_max` (`:165`) and only then falls back. The
`{0,0,0}` sentinel is never returned as a layout.

**The termination guarantee** is `makes_progress` (`:54`): a merge bin may never hold every
sequence the IBF was given, so each level strictly shrinks the input and the recursion depth is
finite by construction. `fallback_layout` (`:64`) backs it up — one user bin per technical bin
while they fit, round robin beyond — which is always a valid layout and always makes progress.

**Defect 4 — merge bin occupancy is capped** in
`code/templates/fast_construct_binning_core.tpp:63`:

```cpp
if (!force && res[b].size() >= bins) return false;
```

A merge bin becomes an IBF with `bins` technical bins one level down, so it must not hold more
user bins than that child can lay out. The `force` path in the fallback still bypasses the cap,
so no sequence is ever dropped.

`generate_hibf` also prints one line per level now, which is enough to see the shape of a run:

```
[hibf] level 0: 1 IBF(s), 128 merge bin(s) to expand
[hibf] level 1: 128 IBF(s), 0 merge bin(s) to expand
```

## Three output format bugs, found while validating

These were not part of the level explosion. They turned up when the layouts were checked
against chopper's parser, and all three are in `code/fast_construct_bin.cpp`.

1. **`:124` — the lower level header was missing its colon.** It wrote
   `fullest_technical_bin_idx17` where `seqan::hibf::layout::layout::read_from` expects
   `fullest_technical_bin_idx:17`. The top level line (`:121`) had the colon; the lower level
   line did not.
2. **`:43` — `@CHOPPER_USER_BINS_END` was written without a trailing newline**, so it ran into
   the `@CHOPPER_CONFIG` that followed it and produced the line
   `@CHOPPER_USER_BINS_END@CHOPPER_CONFIG`.
3. **`:38-44` — the user bin section was emitted in `std::unordered_map` order.**
   `chopper::layout::read_filenames_from` reads that section *positionally*
   (`filenames.emplace_back()`); the `@<id>` on each line is only cross checked under `assert`,
   so in a release build every user bin was silently paired with a file belonging to a different
   one. The section is now emitted by ascending id.

The third one means any layout produced before `a0af981` had scrambled user bin to file
assignments, which is worth knowing before comparing against older benchmark numbers.

## What is still open

* **The two largest datasets keep a small third level** — 2 IBFs for 262144, 17 for 524288.
  The cause is the occupancy cap being exactly `bins`, which leaves a child no slack for a
  sequence that has to be split across two technical bins. For 262144 the geometry is exactly
  tight: every level-1 IBF holds exactly 512 user bins in exactly 512 technical bins
  (512 x 512 = 262144), and 2 of the 512 could not fit. For 524288 the level-1 IBFs hold
  640-768 user bins in 768 bins, and the 16 fullest ran out (one of them twice, hence 17).
  Lowering the cap does not help 262144, because the root has no slack either — fixing that
  last 0.4% / 2.2% means giving `global_bins`
  (`code/templates/fast_construct_generate_hibf.tpp:176`) headroom rather than deriving it from
  `ceil(sqrt(N)/64)*64` alone.
* **1024 and 2048 use only 16 and 32 of their 64 top level technical bins.** One merge bin per
  genome group is the natural clustering for this data, but `binning_core` sorts empty bins to
  the *front*, so the IBF still allocates all 64 with 48 (resp. 32) of them empty — about 4% of
  that layout's technical bins. Sorting empty bins last would drop them, but it breaks the
  `partition_point` logic that derives `split_start` / `merge_start`, so the convention was left
  alone. Note also that 16 is not obviously better than chopper's 64: fewer, fatter top level
  bins may cost query time, so this wants benchmarking rather than assuming.
* **`function_tests` reports 4 failures** — seeding-after-splitting, climbing, splitting, and
  merge/split isolation. These predate the fix: the test output is byte-identical to a build
  from before it (51 TRUE / 4 FALSE both ways).
* **Content line order is not stable across builds.** `write_content` iterates an
  `unordered_map`, so two builds produce semantically identical layouts whose content lines are
  in different order. The parser does not care, but sorting by user bin id would make output
  byte-reproducible and easier to diff.

## Testing

`generate_layout` takes an optional 13th argument, a sketch cache file:

```
./generate_layout <dir> <k> <q> <w> <seed> <fpr> <lvls> <s> <refinements> <hash_funcs> \
    <out> <threads> [sketch_cache]
```

Hashing the FASTA input is by far the slowest step, and it does not change while the layout
algorithm is being worked on. The cache stores `dir_path`, `q`, `k`, `w`, `seed` and `s`
alongside the sketches and refuses to load when any of them differs
(`code/fast_construct_hashing.cpp:97` and `:127`), so a stale cache cannot be picked up
silently. Measured on this data:

| dataset | cold run | from cache |
|---|---:|---:|
| 16384   |  363 s |  13 s |
| 262144  |  983 s |  54 s |

The cached run of 262144 produced a byte-identical layout to the cold one.
