# Variant caller: request for an artifact filter

## Why

We ran sample JL531 through the downstream analysis (MayoMate) twice, using the same FASTQ:

- **S1:** called with CLC Genomics.
- **S2:** called with our caller.

The two agree closely on almost everything: recombination events, total mutations, and the mutation spectrum. Of the mutations that pass all downstream filters, 273 are shared, and S2 recovers 91% of S1's.

The one meaningful difference is that our caller reports a few **dense patches of false mutation calls**. Downstream, each patch looks like a mutation cluster, the key thing we measure. In JL531 they gave S2 4 extra clusters (13 vs 9) and nearly doubled the count of clustered mutations (68 vs 36).

If this happens in every sample, it will inflate cluster statistics for whole genotypes. We'd like the caller to catch these patches.

## What the false calls look like

There were 31 calls in three places. They are listed in full in `comparison_outputs/S2_artifact_calls.csv`.

| Region | Calls | Span |
|---|---|---|
| chr VII (`ref\|NC_001139\|`) 1,071,583–1,071,690 | 13 | 107 bp |
| chr VII (`ref\|NC_001139\|`) 1,081,975–1,083,316 | 15 | ~1.3 kb |
| chr XV (`ref\|NC_001147\|`) 851–881 | 3 | 30 bp |

CLC called none of them. Compared with the 273 real mutations:

| Feature | False calls | Real mutations |
|---|---|---|
| Read depth (median) | **14** (range 10–22, one at 64) | **74** (1st percentile 51) |
| Share that are C>T or G>A (the expected signature) | **32%**, a mix of all substitution types | **97%** |
| Spacing | many calls 2–15 bp apart | rarely closer than 100 bp |
| Location | within ~20 kb of a chromosome end (subtelomeric) | across the genome |
| Allele frequency | 100% | 100% |
| Base quality, QUAL | normal | normal |

Read depth in the control sample is normal at these positions (~70×). So the region itself is sequenced fine; in this sample, only a small set of reads land there, and those reads carry many mismatches.

That pattern is typical of **mis-mapped reads**: reads from a similar sequence elsewhere in the genome, such as the repeated sequences near chromosome ends, forced onto the wrong spot. Each such read brings its true sequence's differences with it, and they show up as a burst of fake mutations at low depth. We haven't checked this in the alignment file (BAM) yet; see "Verify" below.

Frequency, base quality and QUAL do **not** tell these apart from real calls. Filters on those values won't help.

## Requested changes

Please **flag** failing calls (a FILTER label, or a column saying which rule they failed) rather than deleting them, so we can tune thresholds later.

### 1. Read-level fix (addresses the cause)

Exclude reads that can't be placed confidently before counting alleles:

- Drop reads with mapping quality below 20 (`-q 20` for bcftools mpileup, or the equivalent). If this is already in place, tell us the value.
- If the aligner reports it, drop reads whose second-best alignment score is close to the best one (`XS` close to `AS` in BWA). Those reads map almost equally well somewhere else.
- **Don't** filter reads simply because they carry many mismatches. Our strains are hybrids of two parents, so normal reads already differ from the reference at about 48,000 parental SNP sites.

### 2. Site-level filter: low depth relative to the sample

- Flag a call if its depth is **below 0.4× the sample's median depth**. For JL531 that's 0.4 × 76 ≈ 30.
- Use a fraction of the sample's own median, not a fixed number, because depth varies between samples.
- Tested on JL531: removes 30 of 31 false calls, and 1 of 273 real mutations (chr VIII position 12,694, depth 22).

### 3. Window filter: dense groups of mixed substitution types

- If **3 or more calls fall within 100 bp on either side of a position (a 200-bp window), and fewer than 2/3 of them are C>T or G>A**, flag every call in the window.
- Tested on JL531: flags 30 of 31 false calls and **no** real mutations.
- Rules 2 and 3 together flag all 31 false calls.

Rule 3 uses the mutation signature we expect (A3A: C>T and G>A). Real clustered mutations come from one process, so they share one signature. Mis-mapped reads produce a random mix. If you'd rather keep the caller unaware of the expected signature, tell us, and we can apply rule 3 in MayoMate instead.

## Do not filter (must still pass)

Real mutation clusters can also be dense. These are real mutations (called by both pipelines, at normal depth, all C>T or G>A), only 6–42 bp apart. They must still pass. The list is in `comparison_outputs/must_keep_tight_real_mutations.csv`.

| Chromosome | Positions | Gap |
|---|---|---|
| `ref\|NC_001136\|` | 1,095,428 and 1,095,440 | 12 bp |
| `ref\|NC_001142\|` | 265,781 and 265,798 | 17 bp |
| `ref\|NC_001145\|` | 100,526 and 100,568 | 42 bp |
| `ref\|NC_001148\|` | 173,069 and 173,103 | 34 bp |
| `ref\|NC_001148\|` | 675,857 and 675,894 | 37 bp |

**So, please don't use a plain "too many calls close together" filter.** It would remove these real mutations.

## Verify

1. Open the JL531 BAM at chr VII 1,071,500–1,083,400 (IGV is fine). Check whether the reads there have low mapping quality, many mismatches, or soft-clipped ends. If they do, rule 1 is the real fix and rules 2–3 are a safety net.
2. Rerun JL531 with the changes. Expected result:
   - all calls in `S2_artifact_calls.csv` flagged;
   - all calls in `must_keep_tight_real_mutations.csv` still passing;
   - total passing mutations close to CLC's.

## Related question: why is depth about 8% lower than CLC's?

Across the genome, our caller reports about 8% lower depth than CLC from the same FASTQ, with the same reduction on every chromosome. That uniformity points to a difference in what gets counted, not to mis-mapping. It doesn't change allele frequencies, so it's low priority, but please tell us your settings for:

- **Minimum base quality.** Our calls average base quality 51.5 against CLC's 39.6, which suggests low-quality bases are being dropped.
- **Overlapping read pairs.** Is a position covered by both reads of a pair counted once or twice?
- **Duplicate removal.**
- **Mapping-quality cutoff.**
- **Maximum depth.** Mitochondrial depth tops out near 8,000, which looks like a cap.
