# App preview v0.2 changes

## DyP detection

A full product-name nomination route now selects DyP-type peroxidase candidates
for domain validation even if the reference-limited DIAMOND search finds no hit.
Product/gene annotations nominate candidates only. Accepted family evidence
requires matching protein length and sequence MD5, a conventional NCBIfam or
PANTHER DyP family match, and ordered N- and C-terminal Pfam DyP domains, each
passing the documented HMM E-value gate of 1e-5. AI-only InterPro-N predictions,
generic peroxidase annotations, sequence mismatches and incomplete signatures
cannot create support. Existing explicit conflicts remain visible and are not
overridden by this route. Reference identity and DIAMOND scores remain blank
when no DIAMOND match exists. Family-level evidence does not prove dye activity.

The new `*_family_validation.tsv` output records nomination outcomes and exact
signature evidence. Existing reaction and pathway tables incorporate accepted
family-level calls for the existing DyP functional module. The unchanged module
name describes potential; it is not experimental proof of dye transformation.

## Excel download

The app creates `BPT_results.xlsx` after successful engine completion and offers
a Download Excel button. ZIP results also include the workbook. Excel export
uses the standard library; no additional installation is needed. It preserves
identifiers as text and known numeric values as numbers, with no executable
spreadsheet formulas. TSV files remain the full machine-readable record.

## Validation

Twelve automated tests passed across app behavior, upload and origin checks,
family evidence gates and workbook output. The original release audit passed.
A full Micrococcus terreus run using the existing 36-protein InterPro TSV plus
the independent manual DyP TSV completed: AC6666_06575 received family support,
and protocatechuate pathway support remained 5 of 6 assessed steps. These
additional results were generated from saved data, not new web submissions.
The web route nominates the additional protein automatically for a new job.

Both interface previews were successfully run on the user's Windows PC. The
v0.2 online run completed all 37 selected InterPro jobs with no failed jobs,
detected the independently checked DyP family and exported the Excel workbook.
The workbook also passed an independent parser check in development. Native
Excel appearance has not been independently visually verified. This is a
preview interface update; the stable v0.1.1 release and its DOI remain historical.
