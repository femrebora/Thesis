# metabolomics/upstream/ — MetaboAnalystR session records

> **The metabolomics results were removed from the submitted thesis** (unpublished
> data). This directory exists for provenance only.

`metaboanalyst_session_{1,2,3}.R` are saved R history files from the interactive
MetaboAnalystR v4.0.0 sessions that produced the processed tables in
`metabolomics/data/` — normalisation, Welch's t-test, fold-change, PCA/PLS-DA and
the SMPDB/KEGG over-representation analysis.

They are **session transcripts, not runnable scripts**: they contain the calls in
the order they were issued, including exploratory dead ends, and they depend on
MetaboAnalystR's working-directory side effects. Read them to see exactly which
functions and parameters were used; do not expect `Rscript` to execute them
cleanly end to end.

The runnable, cleaned-up figure code is in `metabolomics/scripts/`, which reads
the processed CSVs in `metabolomics/data/` rather than re-deriving them.

## Known issue

`metabolomics/scripts/run_metabolomics_analysis.R` regenerates most figures
correctly but does not terminate on its own — it hangs after the enrichment step,
most likely on a MetaboAnalystR network call. Interrupt it once the figures are
written. Separately, `volcano_enhanced.pdf` regenerates substantially smaller
than the committed version, which suggests content loss in that one figure and
is unresolved. Neither affects the thesis, since this module is excluded from it.
