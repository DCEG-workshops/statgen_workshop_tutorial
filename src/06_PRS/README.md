# Session 6 Polygenic Risk Scores

Compare clumping and thresholding (C+T) with LDpred2 using the same chromosome 22 workshop data in Biowulf OnDemand RStudio. The [R Markdown tutorial](06_PRS_LDpred2.Rmd) contains all student instructions, including package installation, data checks, live fitting, and the prepared backup.

Before class, launch RStudio with R 4.6, 2 CPUs, 8 GB RAM, at least 2 hours, and 5 GB local scratch. In your own clone, open `src/06_PRS/06_PRS_LDpred2.Rmd`, set R's working directory to the source file location, and follow its package installation section. Confirm access to the shared inputs at `/data/DCEG_shared/statgen_workshop_2026/data/workshop6/data/`.

Download the supplied LD matrix once from the [Session 6 data release](https://github.com/DCEG-workshops/statgen_workshop_tutorial/releases/tag/prs-tutorial-2026-10-27). From `src/06_PRS/` in the Terminal:

```bash
curl --fail --location --retry 3 \
  --output 06_PRS_1000G_LD.rds \
  https://github.com/DCEG-workshops/statgen_workshop_tutorial/releases/download/prs-tutorial-2026-10-27/06_PRS_1000G_LD.rds
sha256sum -c SHA256SUMS
```

The 276 MiB LD file is a release asset and is not included by `git clone`. Check the files before editing the Rmd. On macOS, use `shasum -a 256 -c SHA256SUMS`.

| File | Purpose |
|:---|:---|
| `06_PRS_LDpred2.Rmd` | Editable tutorial; defaults to live LDpred2 using supplied LD |
| `06_PRS_LDpred2.prepared.html` | Completed reference report; download and open in a browser |
| `06_PRS_LDpred2_prepared.rds` | Prepared weights and scores for the backup workflow |
| `06_PRS_1000G_LD.rds` | Precomputed LD input, downloaded from the release |
| `SHA256SUMS` | Checksums for the four tutorial files |

Run setup and subsequent chunks in order, or Knit. For the backup, change only `ldpred_mode` to `"prepared"`, save, and follow the Rmd's rerun instructions. Prepared mode needs neither the LD matrix nor `bigsnpr`/`bigsparser`. Leave `rebuild_ld: false` during class; optional LD construction is retained in the instructor appendix.

Raw genotypes/outcomes, installed packages, and generated analysis folders are not included in Git. Results are written to a new subdirectory of `06_PRS_results/` beside your copy of the Rmd.
