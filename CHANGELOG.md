# Changelog

All notable changes to this project are documented in this file.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Every code change gets an entry here and a patch bump (see `CLAUDE.md`).

## [0.0.3] - 2026-09-25

### Added

- Library folders: nested, collapsible folders (state persisted), "New folder", rename,
  delete (contents move up one level), "Move to…" for documents and folders with cycle
  protection. Imports, demo data, mapping results and slices are stored in the current
  folder. Existing libraries are upgraded automatically (`folders` table, `folder_id` column).
- "Extract slice" (selection panel and oligo report): saves the selected columns of all
  sequences, without the reference and with the original names and metadata, as a new
  library document; empty sequences and all-gap columns are dropped
  (`Document::slice_columns`).
- Rename dialogs focus the name field with its text selected.

### Changed

- Identity and conservation graphs never count gaps as identical: they compare against the
  reference (when highlighting disagreements to it) or the column's majority state; a gap or
  missing comparison residue scores 0, and gap/gap pairs no longer count as identical pairs.
  The graphs update when the highlight mode changes.
- Zoomed out, the graph shows the mean of each pixel's columns instead of the minimum.
- Larger folder expand/collapse arrow in the library.

## [0.0.2] - 2026-09-25

Initial implementation of PCR Studio: a desktop sequence viewer for PCR
design and inclusivity evaluation (see `README.md`).

### Added

- Dependencies: `iced` 0.14 (canvas, tokio), `bio` 4, `rusqlite` 0.40 (bundled SQLite),
  `rayon`, `rfd`, `serde`, `serde_json`. Dev profile optimizes dependencies (`opt-level = 3`)
  and this crate lightly (`opt-level = 1`) so debug builds handle large alignments.
- Data model (`model.rs`): rows placed in alignment coordinates with a start column and
  "no coverage" outside it, documents, annotations, parallel per-column statistics,
  IUPAC consensus at selectable thresholds, pairwise identity, conservation, and
  column ↔ reference-position mapping.
- Nucleotide utilities (`seq.rs`): normalization, IUPAC masks, reverse complement,
  GC content, ambiguity expansion.
- FASTA import/export via rust-bio (`fasta.rs`), aligned vs. unaligned detection,
  terminal gaps imported as missing data, `key=value` / `[key=value]` header metadata,
  CSV/TSV metadata tables and name-field splitting.
- Melting temperature (`thermo.rs`): SantaLucia 1998 nearest neighbor with SantaLucia salt
  correction and Mg²⁺/dNTP conversion; matches primer3 `oligotm` to 0.01 °C in tests.
  IUPAC oligos report a Tm range. (The `primer3` crate does not build on Windows/MSVC;
  a pure Rust implementation was chosen instead.)
- Map to reference (`mapper/`): k-mer index with repeat masking, both orientations,
  anchor chaining with rust-bio `sparse::sdpkpp`, gap filling with rust-bio's pairwise
  DP aligner (global on the query, local on the reference), optional clipping of
  overhanging ends, identity/length filters, parallel mapping with progress and cancel,
  contig building that pools insertions into shared columns, and fine tuning that
  re-maps to the consensus until it converges. Unmapped sequences can be saved separately.
- Oligo evaluation (`primer.rs`): Tm, GC, perfect / 1 / 2 / ≥3 mismatch counts, 3' mismatch
  window, no-coverage count, mismatches per position, most common target variants and
  per-group inclusivity; forward or reverse orientation, consensus or reference as source.
- Display list (`display.rs`): sort by name, position, length, identity to reference or
  any metadata field; group by metadata (foldable groups); collapse identical sequences;
  name filter.
- Motif search (`search.rs`): IUPAC aware, 0–3 mismatches, both strands, across gaps, in
  reference, consensus or all sequences.
- SQLite library (`db.rs`) in `pcr_studio.db` next to the executable: documents, rows,
  metadata, annotations and settings.
- Iced user interface (`ui/`): library panel, alignment canvas (overview strip, ruler with
  reference coordinates, consensus, conservation/identity graph, annotation lanes, pinned
  reference, sequence rows), zoom from single bases to the whole genome, adjustable row
  height, scrollbars, wheel/drag/keyboard navigation, jump to position or range,
  highlighting of disagreements to reference or consensus with dots, gap highlighting,
  column selection with the oligo report popup, annotations, map-to-reference, export and
  settings dialogs, background jobs with progress, light and dark themes.
- Synthetic demo data (`demo.rs`): 35 kb reference plus 2,000 related genomes from six clades.
- FASTA files passed on the command line are imported at startup ("Open with").
- `README.md` with usage and architecture notes.

### Changed

- `src/main.rs` starts the application instead of printing "Hello, world!".
- Version bumped to `0.0.2`.

## [0.0.1] - 2026-09-25

### Added

- Initial Rust binary crate `pcr_studio_1` (edition 2024) with `src/main.rs` hello-world entry point.
- `CHANGELOG.md` to document all changes to the codebase.
- `CLAUDE.md` with project overview and working rules (changelog, versioning, git).

### Changed

- Version initialized to `0.0.1` in `Cargo.toml` (was the `cargo new` default `0.1.0`).
