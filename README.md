# DESCARTES Exp2 microbiome

Quarto book for collecting, reviewing and publishing the downstream microbiome
results from `descartes_exp2`.

## Local workflow

Synchronize the latest Markdown reports, figures and tables from the downstream
analysis. The web repository stores figures as `PNG` and tables as `CSV`, even
when the analytical project keeps `SVG/PDF` figures and `TSV/XLSX` tables:

```bash
Rscript scripts/sync_downstream_reports.R
```

Render the website:

```bash
quarto render
```

The rendered site is written to `docs/`, ready for GitHub Pages.

## Source project

Default source:

```text
/mnt/lustre/scratch/nlsas/home/otras/pia/dci/descartes_exp2/custom_downstream_analysis
```

Override it with:

```bash
DESCARTES_DOWNSTREAM=/path/to/custom_downstream_analysis Rscript scripts/sync_downstream_reports.R
```

## Structure

- `_quarto.yml`: Quarto book configuration.
- `chapters/generated/`: synchronized chapters generated from downstream reports.
- `assets/results/`: copied `PNG` figures and `CSV` tables needed by the web book.
- `scripts/sync_downstream_reports.R`: report/assets synchronization script.
- `.github/workflows/pages.yml`: GitHub Pages deployment workflow.
