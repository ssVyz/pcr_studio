# PCR Studio

Desktop application for working with large sequence alignments, with a focus on
PCR design and inclusivity evaluation. Inspired by Geneious Prime's contig viewer
and Map to Reference.

- Rust (edition 2024), UI in [iced](https://iced.rs) 0.14, alignment code from
  [rust-bio](https://rust-bio.github.io), library stored in SQLite.
- Built for tens of kilobases × thousands of sequences: a 35 kb reference with
  2,000 genomes maps in about 4–10 s and stays interactive in the viewer.

## Build and run

```sh
cargo run --release        # recommended for large data sets
cargo run                  # debug build (dependencies are optimized, see Cargo.toml)
cargo test                 # unit tests
cargo test --release -- --ignored bench_genomes --nocapture   # mapping benchmark
```

The document library is the SQLite file `pcr_studio.db` next to the executable
(`target/release/pcr_studio.db` when started with `cargo run --release`).
FASTA files given on the command line are imported at startup.

## Workflow

The menu bar holds all application commands:

| Menu | Entries |
|---|---|
| File | Import FASTA… (Ctrl+O), Import metadata table…, Generate demo data, Export FASTA… (Ctrl+E), Extract selection to library, New folder, Close document (Ctrl+W), Exit |
| Edit | Copy oligo (Ctrl+C), Copy oligo report, Add annotation from selection, Clear selection (Esc), Set selected row as reference, Clear reference, Settings… (Ctrl+,) |
| Tools | Map to reference… (Ctrl+M), Metadata from sequence names…, Oligo report for selection, Find sequence… (Ctrl+F), Go to position… (Ctrl+G) |

Entries that do not apply (e.g. without an open document or selection) are disabled.

1. **Import FASTA** (File menu or Ctrl+O). Aligned FASTA (gapped, or all sequences the
   same length) becomes an alignment; anything else a sequence list. Leading/trailing
   gaps are treated as missing data, not deletions. `key=value` and `[key=value]`
   tokens in FASTA descriptions become metadata. **Demo data** adds a synthetic
   35 kb × 2,000 genome set.
2. **Open** a document from the library (double-click it). Icons distinguish folders,
   sequence lists, alignments and contigs. Organize the library in folders: **New folder**
   creates one inside the selected folder (or next to the selected document); **drag** a
   document or folder onto a folder to move it there, or onto the free space below the
   entries to move it to the top level (hovering a collapsed folder expands it, the list
   scrolls at its edges, Escape cancels, invalid targets are marked red). **Move to…** does
   the same from a list. The chevron or a double-click expands/collapses a folder (remembered). Deleting a folder keeps its contents (they move up one
   level). New documents (imports, mapping results, slices) go to the current folder.
3. **Map to Reference** (Tools menu or Ctrl+M): pick the reference from this or another document,
   choose a sensitivity preset or custom values, optional fine tuning. The result is
   a new contig document with the reference pinned on top.
4. **Evaluate an oligo**: drag across columns (or type a range like `19160-19180` into
   the position box). The oligo report shows 5'→3' sequence, length, Tm, GC, perfect /
   1 / 2 / ≥3 mismatch counts, 3' end mismatches, sequences without coverage,
   mismatches per position, the most common target variants, and, when rows are
   grouped, inclusivity per group. Switch between consensus/reference as the source
   and forward/reverse orientation. Add the selection as an annotation (forward
   primer, reverse primer, probe, region).
5. **Extract slice** (selection panel or oligo report): saves the selected columns of all
   sequences, without the reference and with the original names, as a new library document
   next to the source. Sequences without bases in the range and empty columns are dropped.
6. **Export** (File menu or Ctrl+E) all, shown or selected sequences, optionally only the selected
   columns, gapped or ungapped, with reference and/or consensus.

### Viewer

| Action | How |
|---|---|
| Scroll | Mouse wheel (vertical), Shift+wheel (horizontal), scrollbars, middle-button drag, arrow keys, Page Up/Down, Home/End |
| Zoom | Ctrl+wheel at the cursor, −/+ buttons, Ctrl+=/Ctrl+−, **Fit** (Ctrl+0), **1:1** (readable bases) |
| Row height | *Rows* slider: 2 px rows show hundreds of sequences at once |
| Jump | Position box: reference position (or column when there is no reference); a range selects it |
| Search | IUPAC motif with 0–3 mismatches, both strands, in reference, consensus or all sequences; Enter or ‹ › cycles through hits |
| Overview | Top strip: conservation of the whole alignment; click or drag to navigate |
| Select rows | Click names; Ctrl/Shift for several. Used for "Set selected row as reference", export and mapping subsets |
| Groups | Group by any metadata field; click a group header to fold it |
| Escape | Close dialog/report, then clear the selection |

Display options (right panel): consensus threshold (majority or 50–100 %), graph mode
(**conservation** = share of sequences identical to the comparison residue, on a log scale,
colored ≥99.99 / 99 / 95 / 90 / 75 %; or Geneious-style **pairwise identity**). The comparison
residue is the reference when highlighting disagreements to the reference, otherwise the
column's most common state. Gaps never count as identical: a column that is a gap in the
consensus or reference scores 0, and gap/gap pairs do not count as identical pairs.
Further options: highlighting of disagreements to the reference or the consensus, dots for
identical bases, gap highlighting, base colors, sorting, grouping, collapsing identical sequences and a name
filter. Metadata can be imported from a CSV/TSV table (first column = sequence name) or
split from sequence names at a delimiter.

## Architecture

| Module | Content |
|---|---|
| `model.rs` | Rows (start column + gapped data; outside = no coverage), documents, annotations, column statistics, consensus, identity/conservation, reference coordinates |
| `seq.rs` | Normalization, IUPAC masks, reverse complement, GC, ambiguity expansion |
| `fasta.rs` | FASTA I/O (rust-bio), description metadata, CSV/TSV tables |
| `thermo.rs` | SantaLucia 1998 nearest-neighbor Tm, salt and Mg²⁺ correction (primer3 defaults, validated against primer3) |
| `mapper/` | `index.rs` k-mer index, `align.rs` seed chaining (rust-bio `sdpkpp`) + DP gap filling (rust-bio pairwise), `contig.rs` insertion-aware merge and consensus target, `mod.rs` orchestration and fine tuning |
| `primer.rs` | Oligo/selection evaluation |
| `display.rs` | Sorting, grouping, collapsing, filtering |
| `search.rs` | Motif search |
| `db.rs` | SQLite library (documents, rows, annotations, folders, settings) |
| `demo.rs` | Synthetic data |
| `ui/` | iced application: `mod.rs` state/messages/update, `menu.rs` menu bar and drop-downs, `view.rs` layout and panels, `canvas.rs` alignment renderer, `dialogs.rs`, `state.rs` viewer state, `style.rs`, `job.rs` background jobs |

### Map to reference

Follows the documented behavior of the Geneious mapper (`context/geneious_map_to_reference_context.md`),
simplified to one reference and non-NGS inputs:

1. Index all k-mers of the reference ("index word length"); drop words occurring more
   often than the repeat limit.
2. For each sequence and its reverse complement, collect k-mer hits and chain them with
   rust-bio's sparse dynamic programming (`sdpkpp`); keep the better-scoring orientation.
3. Turn the chain into exact-match segments and fill the gaps between them, and both ends,
   with rust-bio's affine-gap DP aligner (global on the sequence, local on the reference).
   Very large gaps are re-seeded with shorter words. Overhanging ends are clipped
   (optional).
4. Accept alignments above the minimum identity and aligned length.
5. Merge all pairwise alignments into a contig: every reference base is a column,
   insertions at the same position share a block of columns (the reference gets gaps there).
6. Fine tuning (optional, up to N rounds): compute the consensus of the mapped sequences,
   re-map everything to it, and keep the reference aligned to the new target. Stops early
   when the consensus no longer changes. Indels then line up between sequences instead of
   being placed independently per sequence.

Mapping runs in parallel (rayon) and is deterministic.

### Rendering

The viewer is a single iced canvas that draws only the visible window. At readable zoom
each cell is a colored block with a letter. Below one pixel per column, every pixel column
summarizes its columns (the first disagreement wins when highlighting), so a 35 kb × 2,000
alignment can be viewed and scrolled as a whole. Regions are drawn in clipped sub-frames
because iced renders canvas text above all shapes of a frame.
