# Statistical Genetics Workshop — Tutorials

Tutorial materials for the DCEG Statistical Genetics Workshop.

## Workshop years & repository structure

This repository holds the workshop tutorials across multiple years. The
published workshop site and the "Open in Colab" badges link **directly** to
files in this repo at a specific branch, so to keep each past year's links
working we use the following convention:

- **`main`** — the current / upcoming workshop (2026). This is the active
  development line.
- **`YYYY` branches** — frozen snapshots of past workshops. These are never
  rewritten, so links into them stay valid forever.

### Past workshops

| Year | Branch | Materials |
|------|--------|-----------|
| 2023 | [`2023`](https://github.com/DCEG-workshops/statgen_workshop_tutorial/tree/2023) | Google Colab notebooks (`src/*.ipynb`). |

The published 2023 site
(<https://dceg-workshops.github.io/statistical_genetics_workshop/2023/>) links
to the `2023` branch, and the in-notebook "Open in Colab" badges on that branch
resolve to `blob/2023` so they remain pinned to the 2023 materials.

## Current workshop (2026)

Materials for the 2026 workshop are developed on `main`. The tutorials cover:

1. Ancestry & quality control
2. GWAS & meta-analysis
3. Fine-mapping & colocalization
4. Heritability & polygenic risk scores
5. Rare variants
6. Mendelian randomization
7. Multi-ancestry analysis
8. Mosaic chromosomal alterations (mCA)
9. Functional genomics

> **Note:** For 2026 the tutorials are being migrated from Google Colab
> notebooks to R Markdown for use on Biowulf. That conversion is in progress.

## Maintainer note — archiving a workshop year

When a workshop concludes and `main` moves on to the next year:

1. **Freeze the year.** Create a branch from the final commit and push it:
   `git branch YYYY <commit> && git push origin YYYY`.
2. **Make the archive self-consistent.** In that branch, repoint any
   self-referential links (Colab badges, `blob/`/`tree/` links to this repo)
   from `blob/main` to `blob/YYYY`. Leave links to *other* repos untouched.
3. **Repoint the published site.** Update the workshop site's `YYYY` page to
   link to `blob/YYYY` instead of `blob/main` **before** `main` diverges to the
   next year's content.
4. **Record it** in the "Past workshops" table above.
