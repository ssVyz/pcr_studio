# Changelog

All notable changes to this project are documented in this file.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Every code change gets an entry here and a patch bump (see `CLAUDE.md`).

## [0.1.0] - 2026-10-23

### Changed

- Version `0.1.0` (minor bump on request): first minor release, no functional changes since `0.0.17`.

## [0.0.17] - 2026-10-23

### Added

- **Serial annotation…** (Edit menu, and right-click in the viewer): annotates many oligos at once.
  Paste FASTA records; the header becomes the annotation name (records without one become
  “Oligo N”). Each oligo is placed where it has the fewest mismatches on either strand of the
  reference, or of the consensus without a reference (a single sequence is its own consensus).
  Plus-strand matches become forward primers, minus-strand matches reverse primers.
  - Ungapped matching (mismatches only, no indels); matches may span alignment gaps. IUPAC codes
    in the oligo match the bases they include; an N in the template counts as a mismatch, so N
    runs do not attract oligos.
  - **Max. mismatches** (default 3): oligos whose best position has more are listed as not found.
  - Live preview while pasting: position (reference positions, or columns), strand and mismatches
    per oligo. Oligos that fit several positions equally well are flagged; the leftmost (plus
    strand first) is used.
  - Records named like an existing annotation, or like an earlier record, are skipped.
  - Not available for linked (read-only) documents.

## [0.0.16] - 2026-10-22

### Added

- **Linked files**: a new library document kind that points to sequence files on disk instead of
  holding the sequences (for FASTQ files, which are too large for the library). PCR Studio never
  changes a linked file. Own icon (chain), colored by format.
  - **File → Import FASTQ (linked)…** (also the welcome screen, right-click on free library
    space, and FASTQ files given on the command line): pick `.fastq`/`.fq` files (also `.gz`) or
    **Add folder…**. The dialog proposes read sets that can be renamed: R1/R2 files are paired
    (`_R1`/`_R2`, or `_1`/`_2` when exactly these two exist), lanes of one sample form one paired
    set, Nanopore chunk files of one run one set. Options: pair R1/R2 files; one read set per
    sample, per file or pair, or all files together.
  - Linking reads every file once, in parallel (progress, **Cancel**): reads, bases, lengths and
    N50, mean quality, GC and N share; R1 and R2 must have the same reads (count and names) or the
    set is not linked. The technology (Illumina / Nanopore) is detected from the read names, else
    the lengths, and can be changed.
  - **Import FASTA…** now asks whether to **copy the sequences into the library** (default, as
    before) or **link the file**.
  - Every file is remembered by path, size, modification time and a fingerprint (hash of size and
    first and last MiB). Missing (moved, deleted, unreachable) and changed files are found at
    start-up and when a linked document is opened or checked, and marked in the library (broken
    chain, red note).
  - **Linked files…** dialog (right-click, or the button in the right panel): summary, files with
    folder, size and state, **Check files**, **Open in viewer**, **Locate…** for a missing or
    changed file. Locate… accepts only the same file (size and fingerprint) and then also finds
    the other missing files from the same folder, including those of other linked documents.
  - Opening a linked document shows it read-only in the viewer: the first 1,000 reads of a read
    set (pairs: R1 and R2 alternating, metadata *Mate*), a FASTA file completely (up to 100,000
    sequences / 100 Mb). Sorting, the reference, annotations and metadata cannot be saved (menu
    entries disabled, the right panel explains it). Map to reference, Find primers, Extract and
    Export need all sequences: enabled for a completely shown FASTA file, disabled for the first
    reads of a read set.
  - Deleting a linked document removes only the link; the confirmation says so.
- `linked` module (`LinkedSet`, file states, fingerprints, relinking, gzip-aware `open`; `fastq`
  reader that also takes wrapped records and empty reads; `scan` summary pass, R1/R2 check,
  technology detection, preview; `group` file grouping), `linked` table and
  `Db::insert_linked_document` / `linked_set` / `list_linked` / `update_linked`,
  `fasta::read_rows`, `DocKind::Linked`, with tests.
- Dependencies: `flate2` (gzip), `xxhash-rust` (file fingerprints).

### Changed

- FASTA import reads gzip-compressed files, and picking a FASTQ file there explains that FASTQ
  files are linked. The import dialog's filter includes `.gz`.
- Map to reference no longer offers linked documents as the reference (as for text documents).

## [0.0.15] - 2026-10-20

### Added

- **External tools**: PCR Studio can run command-line programs for some features. Without a
  tool everything else works; the features that need it are greyed out (Tools menu) or hidden
  (right-click menu, selection panel).
  - **Edit → Settings → External tools** sets each tool's location: bundled with PCR Studio
    (default; `assets/tools/` next to the executable, or in the source tree for `cargo run`
    and copies of the executable), from `PATH`, or a custom path (**Browse…**). **Test** checks
    it (`--version`, minimum version); saving the settings and starting PCR Studio check all
    tools in the background.
  - Tool runs go on alongside other work (status bar: phase, progress, elapsed time, **Stop**),
    one at a time. No console window opens. Exiting during a run asks first and stops the
    tool, so it does not keep running.
- **Find primers (primersearch)**: Tools menu, viewer right-click menu and **Find primers…** in
  the selection panel (needs primersearch 0.1.1 or newer).
  - Runs primersearch on the selected columns of all, the shown or the selected sequences
    (never the reference); with Collapse gaps on, the columns with gaps in the reference are
    left out, as for Extract slice.
  - The dialog sets mode (no ambiguities, incremental, optimize by mismatch), fixed slice,
    orientation, minimum Tm, conserved 3' bases, the mode's options, primer spacing and
    threads; reaction conditions are prefilled from Settings (oligo nM → µM). It shows how many
    sequences are usable: primersearch drops sequences with a gap, missing coverage or an
    N/ambiguity code in the selection. The last parameters are remembered.
  - The tool's text output is saved as a text document in the library (folder of the source
    document) and shown when the run finishes. A header records document, selection,
    sequences, tool version and path, command line, start/finish time, and maps the output's
    positions (columns passed to primersearch) to reference positions and alignment columns.
    Errors from the tool (e.g. search limits) are shown with its message.
- **Text documents** in the library (new kind, own icon): double-click opens a read-only viewer
  with **Copy all** and **Save as .txt…**. They can be renamed, moved and deleted like other
  documents and are not offered as mapping references.
- `tools` module (`ToolSpec`, `ToolConfig`, `check`, `run::run`, `primersearch`),
  `Document::slice_rows` (slice of chosen rows with the alignment column of every slice column),
  `Db::insert_text_document` / `text_content` (`texts` table), with tests; an ignored test runs
  the bundled primersearch (`cargo test -- --ignored bundled_primersearch`).

### Changed

- Closing the window or File → Exit while a tool runs asks for confirmation.
- `.gitignore`: `assets/` (bundled tool binaries are not committed).

## [0.0.14] - 2026-10-13

### Added

- **Sorting sequences** from the viewer's right-click menu and the Edit menu: by name, by
  length (longest first), by the base in the right-clicked column, and by the bases in the
  selected columns (only the shown ones while gaps are collapsed). Bases sort as A, C, G, T,
  other IUPAC codes, gap, no coverage. Sequence lists and alignments/contigs alike.
  - Always sorts all sequences of the document, including those hidden by the name filter
    or folded groups; the reference stays first.
  - The order is saved in the library: reopening, Export FASTA (all sequences) and Extract
    slice keep it. The order from the import is not kept separately.
  - The menu labels show the direction (e.g. "Sort by name (A→Z)"); choosing the same sort
    again right away reverses it. The sort is stable, so ties keep their previous order and
    sorts can be combined.
  - Selected rows and search hits stay with their sequences. A view sort set in the panel is
    reset to "Document order" so the new order is visible (the status bar says so).
- `display::RowSort`, `display::sort_order`, `display::base_rank`, `Db::reorder_rows` (moves
  rows cycle by cycle, about 0.25 s for 2,000 × 35 kb), `OpenDoc::apply_order`, with tests.

### Changed

- Name and metadata comparisons compare runs of digits by value ("seq_9" before "seq_10"),
  also in the panel's Sort and the group order.
- The panel's first Sort option is called "Document order" (was "Original order"), since
  the stored order can now change.
- The viewer's right-click menu is slightly wider to fit the sort entries.

## [0.0.13] - 2026-10-13

### Fixed

- The tooltips of **Collapse gaps** and of **Extract slice** (while gaps are collapsed) showed
  a wide run of spaces in the middle of the sentence.
- The "No metadata yet" hint in the right panel pointed to the toolbar; it now names the File
  and Tools menus.

### Changed

- README: no longer refers to Geneious Prime; describes the mapper on its own terms and the
  pairwise identity graph as the mean identity over all pairs of sequences in a column.
  Added the oligo report's N/ambiguity count, report copying, the Settings dialog, the
  display toggles and `ui/icons.rs`.

## [0.0.12] - 2026-10-13

### Added

- **Collapse gaps** (Display section of the right panel, off by default, saved in the
  settings): hides the columns where the reference has a gap or no coverage, i.e. insertions
  in the other sequences, so the alignment is shown in reference coordinates. Needs a
  reference; without one the option is disabled and says so. The panel shows how many columns
  are hidden, and a ▼ on the ruler marks each place where columns are hidden. Toggling keeps
  the column at the center of the view in place, and a view zoomed to fit stays fitted.
- **Extract slice** leaves out the columns with gaps in the reference while gaps are
  collapsed (the new document's description and the status message say so).
- Export FASTA: **Leave out columns with gaps in the reference** (needs a reference; preset
  while gaps are collapsed), for the whole alignment or the selected columns.
- The selection summary in the right panel shows how many of the selected columns are hidden.
- `model::ColumnMap` (display ↔ alignment columns), `model::pick_columns` and a column
  filter for `Document::slice_columns`, with tests.

### Changed

- The viewer's zoom and scroll position count displayed columns. Selections, annotations,
  search hits and positions still refer to alignment columns; dragging an annotation edge or
  selecting while gaps are collapsed snaps to the shown columns.
- The oligo report still evaluates every column of a selection, including hidden ones.

## [0.0.11] - 2026-10-13

### Added

- **Right-click menu in the sequence viewer**: Add annotation…, Oligo report, Copy oligo,
  Extract selection to library, Clear selection. On an annotation it also offers Annotation
  properties… and Delete annotation…; the right-click selects that annotation.
- **Annotations can be selected**: clicking one (in the viewer or in the side panel list)
  highlights it with grips on both edges and makes its range the column selection. Primer
  annotations set the report orientation, as before.
- **Resizing**: drag either edge of an annotation in the viewer, snapping to columns, with
  a minimum length of one column. The cursor shows a resize arrow over an edge. The new
  range is saved when the mouse is released.
- **Annotation properties dialog** (double-click, right-click menu, Edit menu, ✎ in the
  side panel): name, type and range. The range is entered as 1-based reference positions
  when the document has a reference, otherwise as alignment columns; the dialog shows the
  resulting columns and rejects invalid ranges. A field left unchanged keeps its exact
  column, so annotations that start or end in insertion columns are not shifted. New
  annotations open the same dialog pre-filled from the selection (Add).
- **Deleting an annotation** (right-click menu, Edit menu, Delete key, ✕ in the side panel,
  Delete… in the dialog) asks for confirmation.
- **"Open automatically"** checkbox next to the Oligo report button in the display options
  panel (saved in the settings, on by default). When it is off, selecting columns, jumping
  to a range or clicking an annotation no longer opens the report.
- Edit menu: Add annotation from selection…, Annotation properties…, Delete annotation… (Del).
- `Db::update_annotation`, `RefCoords::range` / `RefCoords::columns`, with tests.

### Changed

- The oligo report no longer contains the annotation name, type and "Add annotation"
  controls.
- An annotation's stored oligo (note) is updated when its range or type changes.

## [0.0.10] - 2026-10-13

### Added

- **New sequence…** (File menu, Ctrl+N, right-click on free library space or on a folder,
  welcome screen): a dialog with name, description and a multi-line sequence box for typed
  or pasted text. Accepts a plain sequence or FASTA. Spaces, line breaks and position numbers
  (GenBank/EMBL listings) are ignored, U becomes T, and characters that are not nucleotide
  codes become N, as on FASTA import. The dialog says how many characters were converted.
  A single FASTA header fills in name and description; several records create one sequence
  list, or an alignment when gapped or all the same length, like a FASTA import. The dialog
  shows length, GC and ambiguous bases, or the records found. The new document goes to the
  right-clicked or current folder and opens.
- **Edit sequence…** (right-click on a document, Edit menu) for documents with a single
  sequence; disabled for other documents. Edits the name, description and bases in the same
  dialog. Renaming there renames the document and its sequence. Metadata is kept.
  Annotations keep their positions; any that end past the new sequence end are removed, and
  the dialog warns before saving. An open document is reloaded in place.
- `fasta::parse_pasted`, `fasta::document_from_pasted`, `fasta::pasted_kind` and
  `Db::replace_single_row`, with tests.

### Changed

- FASTA import and pasted text share one rule for detecting alignments
  (`fasta::alignment_kind`).

## [0.0.9] - 2026-09-28

### Added

- Library multi-selection: Ctrl+click adds or removes entries; dragging any selected entry
  moves the whole selection (the floating label shows "N items"), and **Move to…** moves all
  selected entries too. Entries whose folder is also selected move with that folder; the drop
  is refused if a selected folder would go into itself. A plain click inside the selection
  keeps it for dragging and selects just that entry on release; clicking free space clears
  the selection.
- Right-click on free space in the library: New folder (top level) and Import FASTA….
  Folders keep "New subfolder" in their right-click menu.

### Changed

- The app tracks keyboard modifier state (for Ctrl+click in the library).

## [0.0.8] - 2026-09-28

### Added

- Resizable panels: drag the edges between the library, the sequence view and the display
  options panel (library 200–600 px, options 260–520 px; the sequence view always keeps at
  least 420 px and panels shrink if the window gets too narrow). Double-click an edge to
  reset that panel. The right edge of the sequence-name column in the alignment view can be
  dragged too (120–600 px). All widths are saved in the settings.
- **View** menu: "Display options panel" (check mark, Ctrl+Shift+D) hides/shows the panel on
  the right — showing it again restores its default width — and "Reset panel sizes".
- "Hide" button at the top of the display options panel.
- Menu entries can carry check marks (`menu::Entry::check`).

## [0.0.7] - 2026-09-28

### Changed

- The main window starts maximized (restoring it returns to 1500 × 920).

## [0.0.6] - 2026-09-28

### Added

- Right-click menu on library entries (selects the entry, opens at the cursor and stays
  inside the window; click outside or Escape closes it): documents — Open, Rename…, Delete…;
  folders — Expand/Collapse, New subfolder, Rename…, Delete folder…. Destructive entries are
  shown in red (`menu::Entry::danger`).
- Delete confirmations show what is deleted (type, sequences, length; folder contents that
  are kept), a highlighted "This cannot be undone" notice, and whether an open document will
  be closed. After deleting, the status bar reports what was deleted.
- The window size is tracked (for placing context menus).

## [0.0.5] - 2026-09-28

### Added

- Library icons (inline SVG, `ui/icons.rs`; enables iced's `svg` feature): open/closed
  folders, sequence lists, alignments and contigs, plus SVG expand/collapse chevrons.
- Drag and drop in the library: drag documents or folders onto a folder to move them, or
  onto the free space below the entries to move them to the top level. A floating label shows
  the item and the destination; the target folder is outlined (red and "cannot move here" for
  invalid moves such as a folder into its own subfolder). Hovering a collapsed folder while
  dragging expands it, the list auto-scrolls at its edges, Escape cancels. Drags start after
  a few pixels of movement, so clicking still selects.

### Changed

- Library entries: double-click opens a document or expands/collapses a folder (replaces
  "click the selected entry again"); folder names are bold; hover highlighting.
- Library panel widened from 260 to 290 px.

## [0.0.4] - 2026-09-28

### Added

- Menu bar with drop-down menus (`ui/menu.rs`) replacing the flat toolbar buttons:
  - **File**: Import FASTA…, Import metadata table…, Generate demo data, Export FASTA…,
    Extract selection to library, New folder, Close document, Exit.
  - **Edit**: Copy oligo, Copy oligo report, Add annotation from selection, Clear selection,
    Set selected row as reference, Clear reference, Settings….
  - **Tools**: Map to reference…, Metadata from sequence names…, Oligo report for selection,
    Find sequence…, Go to position….
  Entries show their shortcuts and are disabled when they do not apply; clicking outside or
  Escape closes a menu, hovering another title while one is open switches to it.
- Commands and shortcuts: Close document (Ctrl+W), Find sequence (Ctrl+F, focuses the search
  field), Go to position (Ctrl+G, focuses the position field), Settings (Ctrl+,), Exit.
- The name of the open document is shown at the right of the menu bar.
- Expanded empty library folders show "(empty)".

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
