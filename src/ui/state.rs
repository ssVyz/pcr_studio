//! State of the open document and the viewer.

use crate::display::{self, DisplayOptions, Item, RowSort};
use crate::linked::LinkedSet;
use crate::model::{Annotation, ColumnMap, ColumnStats, ConsensusThreshold, DocInfo, Document, RefCoords};
use crate::primer::{OligoSource, PrimerOptions, PrimerReport};
use crate::search::{Hit, Scope};
use crate::seq::GAP;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::sync::Arc;

/// What to compare bases against when highlighting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Highlight {
    None,
    Reference,
    Consensus,
}

impl std::fmt::Display for Highlight {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Highlight::None => "No highlighting",
            Highlight::Reference => "Disagreements to reference",
            Highlight::Consensus => "Disagreements to consensus",
        })
    }
}

/// What the graph under the consensus shows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GraphMode {
    /// Share of sequences matching the most common residue, log-scaled.
    Conservation,
    /// Mean pairwise identity (Geneious style).
    Identity,
}

impl std::fmt::Display for GraphMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            GraphMode::Conservation => "Conservation (% matching consensus)",
            GraphMode::Identity => "Pairwise identity",
        })
    }
}

/// How reads with base qualities (FASTQ) are shown.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReadDisplay {
    #[default]
    Sequence,
    /// A sequence lane and a quality lane under it.
    SequenceQuality,
    Quality,
}

impl ReadDisplay {
    pub const ALL: [ReadDisplay; 3] = [ReadDisplay::Sequence, ReadDisplay::SequenceQuality, ReadDisplay::Quality];
}

impl std::fmt::Display for ReadDisplay {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            ReadDisplay::Sequence => "Sequence",
            ReadDisplay::SequenceQuality => "Sequence + quality",
            ReadDisplay::Quality => "Quality",
        })
    }
}

/// Persisted display preferences.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ViewPrefs {
    pub highlight: Highlight,
    pub use_dots: bool,
    pub highlight_gaps: bool,
    pub show_consensus: bool,
    pub show_identity: bool,
    pub show_annotations: bool,
    pub show_overview: bool,
    /// Width of the sequence-name column in the alignment view.
    pub name_w: f32,
    pub threshold: ConsensusThreshold,
    pub color_bases: bool,
    pub graph: GraphMode,
    /// Open the oligo report whenever columns are selected.
    pub auto_report: bool,
    /// Hide the columns where the reference has a gap (or no coverage).
    pub collapse_gaps: bool,
    /// Sequence and/or base qualities of FASTQ reads.
    pub read_display: ReadDisplay,
}

impl Default for ViewPrefs {
    fn default() -> Self {
        ViewPrefs {
            highlight: Highlight::Reference,
            use_dots: true,
            highlight_gaps: true,
            show_consensus: true,
            show_identity: true,
            show_annotations: true,
            show_overview: true,
            name_w: NAME_W_DEFAULT,
            threshold: ConsensusThreshold::Majority,
            color_bases: true,
            graph: GraphMode::Conservation,
            auto_report: true,
            collapse_gaps: false,
            read_display: ReadDisplay::Sequence,
        }
    }
}

pub const NAME_W_DEFAULT: f32 = 250.0;
pub const NAME_W_MIN: f32 = 120.0;
pub const NAME_W_MAX: f32 = 600.0;

/// Zoom level where letters are comfortable; shown as 100%.
pub const BASE_COL_W: f32 = 12.0;
/// Space under each read shown with sequence and quality lanes.
const READ_GAP: f32 = 3.0;
pub const MAX_COL_W: f32 = 32.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Viewport {
    /// Pixels per alignment column.
    pub col_w: f32,
    /// Pixels per sequence row.
    pub row_h: f32,
    /// First visible column (fractional).
    pub scroll_x: f32,
    /// Vertical scroll offset of the rows area in pixels.
    pub scroll_y: f32,
    /// Size of the sequence area (without names, headers, scrollbars).
    pub seq_w: f32,
    pub rows_h: f32,
    /// Full canvas size last reported by the canvas.
    pub canvas: (f32, f32),
}

impl Default for Viewport {
    fn default() -> Self {
        Viewport { col_w: BASE_COL_W, row_h: 16.0, scroll_x: 0.0, scroll_y: 0.0, seq_w: 800.0, rows_h: 400.0, canvas: (0.0, 0.0) }
    }
}

impl Viewport {
    pub fn visible_cols(&self) -> f32 {
        self.seq_w / self.col_w
    }
    pub fn min_col_w(&self, width: usize) -> f32 {
        (self.seq_w / width.max(1) as f32).min(BASE_COL_W).max(0.0005)
    }
    /// Clamps zoom and scroll; `content_h` is the height of all rows in pixels.
    pub fn clamp(&mut self, width: usize, content_h: f32) {
        self.col_w = self.col_w.clamp(self.min_col_w(width), MAX_COL_W);
        let max_x = (width as f32 - self.visible_cols()).max(0.0);
        self.scroll_x = self.scroll_x.clamp(0.0, max_x);
        let max_y = (content_h - self.rows_h + self.row_h * 0.5).max(0.0);
        self.scroll_y = self.scroll_y.clamp(0.0, max_y);
    }
    /// Scrolls horizontally so that column `c` is visible (centered if needed).
    pub fn reveal_col(&mut self, c0: usize, c1: usize) {
        let vis = self.visible_cols();
        if (c0 as f32) < self.scroll_x || (c1 as f32) > self.scroll_x + vis {
            let mid = (c0 + c1) as f32 / 2.0;
            self.scroll_x = mid - vis / 2.0;
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Selection {
    pub anchor: usize,
    pub c0: usize,
    /// Exclusive.
    pub c1: usize,
}

impl Selection {
    pub fn new(a: usize, b: usize) -> Selection {
        Selection { anchor: a, c0: a.min(b), c1: a.max(b) + 1 }
    }
    pub fn len(&self) -> usize {
        self.c1 - self.c0
    }
}

pub struct SearchState {
    pub query: String,
    pub scope: Scope,
    pub mismatches: u8,
    pub hits: Vec<Hit>,
    pub current: Option<usize>,
    pub searched: bool,
}

impl Default for SearchState {
    fn default() -> Self {
        SearchState { query: String::new(), scope: Scope::Consensus, mismatches: 0, hits: Vec::new(), current: None, searched: false }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hover {
    pub col: usize,
    pub item: Option<usize>,
    /// Hovered pinned reference row.
    pub on_reference: bool,
}

pub struct OpenDoc {
    pub info: DocInfo,
    pub doc: Arc<Document>,
    pub reference: Option<usize>,
    pub annotations: Vec<Annotation>,
    pub stats: ColumnStats,
    pub consensus: Vec<u8>,
    pub identity: Vec<Option<f32>>,
    pub conservation: Vec<Option<f32>>,
    pub ref_coords: Option<RefCoords>,
    /// Columns shown in the viewer; the viewport counts display columns.
    pub cols: ColumnMap,
    pub row_identity: Option<Vec<f32>>,
    pub meta_keys: Vec<String>,
    pub display: DisplayOptions,
    pub items: Vec<Item>,
    /// Number of group header items before each item (len = items + 1).
    pub group_prefix: Vec<u32>,
    pub vp: Viewport,
    pub selection: Option<Selection>,
    pub selected_rows: HashSet<usize>,
    pub last_clicked_item: Option<usize>,
    pub report: Option<PrimerReport>,
    pub popup_open: bool,
    pub search: SearchState,
    pub jump_text: String,
    pub hover: Option<Hover>,
    /// Selected annotations (ids); the last one is the current one, whose
    /// columns are selected (see [`OpenDoc::selected_annotation`]).
    pub selected_annotations: Vec<i64>,
    /// Where a Shift+click range of annotations starts.
    pub annotation_anchor: Option<i64>,
    /// Type of the last annotation added (default for the next one).
    pub annotation_kind: crate::model::AnnotationKind,
    /// Last permanent sort and whether it was descending (repeating it reverses).
    pub last_sort: Option<(RowSort, bool)>,
    /// Set for linked files: the document is read-only (nothing is saved).
    pub linked: Option<LinkedView>,
    /// How reads are shown; always `Sequence` without base qualities.
    pub read_display: ReadDisplay,
}

/// A linked document in the viewer.
pub struct LinkedView {
    pub set: LinkedSet,
    /// All records of the files are shown, not only the first ones.
    pub complete: bool,
    /// The reads shown are a random sample (otherwise the first ones).
    pub random: bool,
    /// Phred qualities of the bases of each row (FASTQ; empty for FASTA).
    pub qual: Vec<Vec<u8>>,
}

impl OpenDoc {
    pub fn new(info: DocInfo, doc: Document, annotations: Vec<Annotation>, prefs: &ViewPrefs) -> OpenDoc {
        let reference = info.reference.filter(|&r| r < doc.rows.len());
        let meta_keys = doc.meta_keys();
        let mut od = OpenDoc {
            info,
            doc: Arc::new(doc),
            reference,
            annotations,
            stats: ColumnStats::default(),
            consensus: Vec::new(),
            identity: Vec::new(),
            conservation: Vec::new(),
            ref_coords: None,
            cols: ColumnMap::default(),
            row_identity: None,
            meta_keys,
            display: DisplayOptions::default(),
            items: Vec::new(),
            group_prefix: vec![0],
            vp: Viewport::default(),
            selection: None,
            selected_rows: HashSet::new(),
            last_clicked_item: None,
            report: None,
            popup_open: false,
            search: SearchState::default(),
            jump_text: String::new(),
            hover: None,
            selected_annotations: Vec::new(),
            annotation_anchor: None,
            annotation_kind: crate::model::AnnotationKind::ForwardPrimer,
            last_sort: None,
            linked: None,
            read_display: ReadDisplay::Sequence,
        };
        od.recompute_stats(prefs);
        od.rebuild_items();
        if od.reference.is_some() {
            od.search.scope = Scope::Reference;
        }
        od
    }

    pub fn width(&self) -> usize {
        self.doc.width
    }

    /// Direction of the next sort by `key`: the opposite of the last one when
    /// it is repeated, otherwise the key's default.
    pub fn next_sort_descending(&self, key: &RowSort) -> bool {
        match &self.last_sort {
            Some((k, d)) if k == key => !d,
            _ => key.default_descending(),
        }
    }

    /// Reorders the sequences: row `order[i]` becomes row `i`. Row indices held
    /// by the viewer (reference, selected rows, search hits) follow their rows.
    pub fn apply_order(&mut self, order: &[usize]) {
        let mut new_of = vec![0; order.len()];
        for (new, &old) in order.iter().enumerate() {
            new_of[old] = new;
        }
        let doc = Arc::make_mut(&mut self.doc);
        let mut rows: Vec<Option<crate::model::Row>> = std::mem::take(&mut doc.rows).into_iter().map(Some).collect();
        doc.rows = order.iter().map(|&o| rows[o].take().expect("order is a permutation")).collect();
        self.reference = self.reference.map(|r| new_of[r]);
        self.info.reference = self.reference;
        self.selected_rows = self.selected_rows.iter().map(|&r| new_of[r]).collect();
        if let Some(ids) = &mut self.row_identity {
            *ids = order.iter().map(|&o| ids[o]).collect();
        }
        let remap = |h: Hit| Hit { row: h.row.map(|r| new_of[r]), ..h };
        let current = self.search.current.and_then(|i| self.search.hits.get(i).copied()).map(remap);
        for h in &mut self.search.hits {
            *h = remap(*h);
        }
        self.search.hits.sort_by_key(|h| (h.c0, h.row, h.reverse));
        self.search.current = current.and_then(|c| self.search.hits.iter().position(|h| *h == c));
        self.last_clicked_item = None;
        self.hover = None;
        self.rebuild_items();
    }

    /// Number of columns in the viewer (fewer than `width` with gaps collapsed).
    pub fn view_width(&self) -> usize {
        self.cols.len()
    }

    /// Flags of the columns slices and exports keep: the reference's bases
    /// while gaps are collapsed, otherwise `None` (all columns).
    pub fn kept_columns(&self) -> Option<&[bool]> {
        self.ref_coords.as_ref().filter(|_| self.cols.is_subset()).map(|rc| rc.is_base.as_slice())
    }

    /// The shown columns of `[c0, c1)` (all of them when none is shown).
    pub fn shown_columns(&self, c0: usize, c1: usize) -> Vec<usize> {
        let shown: Vec<usize> = (self.cols.disp(c0)..self.cols.disp(c1)).map(|d| self.cols.col(d)).collect();
        if shown.is_empty() { (c0..c1).collect() } else { shown }
    }

    /// Shows all columns, or (`collapse` with a reference) only those where
    /// the reference has a base. The column at the center of the view stays
    /// there; a view zoomed to fit the width stays fitted.
    pub fn update_columns(&mut self, collapse: bool) {
        let center = (self.vp.scroll_x + self.vp.visible_cols() / 2.0).max(0.0);
        let c = self.cols.col(center as usize);
        let fit = !self.cols.is_empty() && self.vp.col_w <= self.vp.min_col_w(self.cols.len()) * 1.001;
        self.cols = match (&self.ref_coords, collapse) {
            (Some(rc), true) => ColumnMap::only(&rc.is_base),
            _ => ColumnMap::all(self.doc.width),
        };
        if fit {
            self.vp.col_w = self.vp.min_col_w(self.cols.len());
        }
        self.vp.scroll_x = self.cols.disp(c) as f32 + center.fract() - self.vp.visible_cols() / 2.0;
        self.clamp_view();
    }

    /// Scrolls horizontally so that columns `[c0, c1)` are visible (centered if needed).
    pub fn reveal_cols(&mut self, c0: usize, c1: usize) {
        let d0 = self.cols.disp(c0);
        self.vp.reveal_col(d0, self.cols.disp(c1).max(d0 + 1));
        self.clamp_view();
    }

    /// Centers the view on columns `[c0, c1)`.
    pub fn center_on(&mut self, c0: usize, c1: usize) {
        let mid = (self.cols.disp(c0) + self.cols.disp(c1)) as f32 / 2.0;
        self.vp.scroll_x = mid - self.vp.visible_cols() / 2.0;
        self.clamp_view();
    }

    /// The per-column values shown in the graph for the chosen mode.
    pub fn graph_track(&self, mode: GraphMode) -> &[Option<f32>] {
        match mode {
            GraphMode::Conservation => &self.conservation,
            GraphMode::Identity => &self.identity,
        }
    }

    /// Recomputes column statistics, consensus, identity and reference coordinates.
    pub fn recompute_stats(&mut self, prefs: &ViewPrefs) {
        self.stats = ColumnStats::compute(&self.doc, self.reference);
        self.consensus = self.stats.consensus(prefs.threshold, b' ');
        self.ref_coords = self.reference.map(|r| RefCoords::new(&self.doc.rows[r], self.doc.width));
        self.update_columns(prefs.collapse_gaps);
        self.recompute_graphs(prefs);
        self.row_identity = None;
    }

    /// Identity/conservation tracks against the comparison sequence: the
    /// reference when highlighting disagreements to it, otherwise the
    /// per-column majority state. Gaps never count as identical.
    pub fn recompute_graphs(&mut self, prefs: &ViewPrefs) {
        let stats = &self.stats;
        let (identity, conservation) = match (prefs.highlight, self.reference) {
            (Highlight::Reference, Some(r)) => {
                let row = &self.doc.rows[r];
                stats.graph_tracks(|c| row.at(c))
            }
            _ => stats.graph_tracks(|c| stats.majority_state(c)),
        };
        self.identity = identity;
        self.conservation = conservation;
    }

    pub fn recompute_consensus(&mut self, prefs: &ViewPrefs) {
        self.consensus = self.stats.consensus(prefs.threshold, b' ');
    }

    /// Per-row identity to the reference over the columns both cover.
    pub fn ensure_row_identity(&mut self) {
        if self.row_identity.is_some() {
            return;
        }
        let Some(r) = self.reference else {
            self.row_identity = Some(vec![0.0; self.doc.rows.len()]);
            return;
        };
        let doc = &self.doc;
        let rref = &doc.rows[r];
        let ids: Vec<f32> = doc
            .rows
            .par_iter()
            .map(|row| {
                let (a, b) = (row.start.max(rref.start), row.end().min(rref.end()));
                let (mut same, mut total) = (0u32, 0u32);
                for c in a..b {
                    let (x, y) = (row.data[c - row.start], rref.data[c - rref.start]);
                    if x == GAP && y == GAP {
                        continue;
                    }
                    total += 1;
                    same += (x == y) as u32;
                }
                if total == 0 { 0.0 } else { same as f32 / total as f32 }
            })
            .collect();
        self.row_identity = Some(ids);
    }

    pub fn rebuild_items(&mut self) {
        if self.display.sort == display::SortKey::IdentityToReference {
            self.ensure_row_identity();
        }
        self.items = display::build(&self.doc, self.reference, &self.display, self.row_identity.as_deref());
        self.group_prefix = Vec::with_capacity(self.items.len() + 1);
        let mut g = 0u32;
        self.group_prefix.push(0);
        for it in &self.items {
            g += matches!(it, Item::Group { .. }) as u32;
            self.group_prefix.push(g);
        }
        self.clamp_view();
    }

    /// Height of group header rows: readable even when sequence rows are compact.
    pub fn group_h(&self) -> f32 {
        self.vp.row_h.max(18.0)
    }

    /// Height of a sequence item: one lane, or a sequence and a quality lane.
    pub fn seq_h(&self) -> f32 {
        match self.read_display {
            ReadDisplay::SequenceQuality => 2.0 * self.vp.row_h + READ_GAP,
            _ => self.vp.row_h,
        }
    }

    /// Top of display item `i` in content pixels (`i == items.len()` gives the total height).
    pub fn item_top(&self, i: usize) -> f32 {
        let g = self.group_prefix[i.min(self.items.len())] as f32;
        (i as f32 - g) * self.seq_h() + g * self.group_h()
    }

    pub fn item_h(&self, i: usize) -> f32 {
        match self.items.get(i) {
            Some(Item::Group { .. }) => self.group_h(),
            _ => self.seq_h(),
        }
    }

    pub fn content_h(&self) -> f32 {
        self.item_top(self.items.len())
    }

    /// Display item at content pixel `y`.
    pub fn item_at_y(&self, y: f32) -> Option<usize> {
        if y < 0.0 || y >= self.content_h() {
            return None;
        }
        let (mut lo, mut hi) = (0usize, self.items.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.item_top(mid + 1) <= y {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        (lo < self.items.len()).then_some(lo)
    }

    pub fn clamp_view(&mut self) {
        let h = self.content_h();
        self.vp.clamp(self.cols.len(), h);
    }

    /// Scrolls vertically so that item `i` is visible.
    pub fn reveal_item(&mut self, i: usize) {
        let (top, h) = (self.item_top(i), self.item_h(i));
        if top < self.vp.scroll_y || top + h > self.vp.scroll_y + self.vp.rows_h {
            self.vp.scroll_y = top - self.vp.rows_h / 2.0;
        }
        self.clamp_view();
    }

    /// Changes the row height, keeping the top visible item in place.
    pub fn set_row_h(&mut self, h: f32) {
        self.keep_top(|od| od.vp.row_h = h);
    }

    /// Shows reads as sequence, quality or both, keeping the top visible item
    /// in place; documents without base qualities always show sequence.
    pub fn set_read_display(&mut self, mode: ReadDisplay) {
        let mode = if self.qualities().is_some() { mode } else { ReadDisplay::Sequence };
        self.keep_top(|od| od.read_display = mode);
    }

    /// Applies a change of item heights, keeping the top visible item in place.
    fn keep_top(&mut self, change: impl FnOnce(&mut OpenDoc)) {
        let top = self.item_at_y(self.vp.scroll_y);
        change(self);
        if let Some(i) = top {
            self.vp.scroll_y = self.item_top(i);
        }
        self.clamp_view();
    }

    /// Phred qualities of the bases of every row (linked FASTQ reads).
    pub fn qualities(&self) -> Option<&[Vec<u8>]> {
        self.linked.as_ref().map(|l| l.qual.as_slice()).filter(|q| !q.is_empty())
    }

    /// Phred quality of row `row` at alignment column `col`.
    pub fn quality_at(&self, row: usize, col: usize) -> Option<u8> {
        let q = self.qualities()?.get(row)?;
        q.get(col.checked_sub(self.doc.rows.get(row)?.start)?).copied()
    }

    /// Takes over the display options and the horizontal zoom and position of
    /// `old`, another set of reads of the same document.
    pub fn keep_view(&mut self, old: &OpenDoc) {
        let mut display = old.display.clone();
        if display.group_by.as_ref().is_some_and(|k| !self.meta_keys.contains(k)) {
            display.group_by = None;
        }
        if matches!(&display.sort, display::SortKey::Meta(k) if !self.meta_keys.contains(k)) {
            display.sort = display::SortKey::Original;
        }
        self.display = display;
        self.vp.col_w = old.vp.col_w;
        self.vp.scroll_x = old.vp.scroll_x;
        self.rebuild_items();
    }

    /// Full-width template (consensus or reference with blanks as gaps).
    pub fn template_row(&self, source: OligoSource) -> Vec<u8> {
        match (source, self.reference) {
            (OligoSource::Reference, Some(r)) => {
                let row = &self.doc.rows[r];
                (0..self.doc.width).map(|c| row.at(c).unwrap_or(GAP)).collect()
            }
            _ => self.consensus.iter().map(|&b| if b == b' ' { GAP } else { b }).collect(),
        }
    }

    pub fn evaluate_selection(&mut self, opts: &PrimerOptions) {
        let Some(sel) = self.selection else {
            self.report = None;
            return;
        };
        let template = self.template_row(opts.source);
        self.report = Some(crate::primer::evaluate(
            &self.doc,
            self.reference,
            &template,
            sel.c0,
            sel.c1,
            opts,
            self.display.group_by.as_deref(),
        ));
    }

    /// 1-based reference range covered by columns `[c0, c1)`, if there is a reference.
    pub fn ref_range(&self, c0: usize, c1: usize) -> Option<(u32, u32)> {
        self.ref_coords.as_ref()?.range(c0, c1)
    }

    pub fn annotation(&self, id: i64) -> Option<&Annotation> {
        self.annotations.iter().find(|a| a.id == id)
    }

    /// The current annotation: the last one selected.
    pub fn selected_annotation(&self) -> Option<i64> {
        self.selected_annotations.last().copied()
    }

    pub fn clear_annotation_selection(&mut self) {
        self.selected_annotations.clear();
        self.annotation_anchor = None;
    }

    /// Updates the annotation selection for a click on annotation `id`: alone,
    /// Ctrl adds or removes it, Shift selects from the anchor to it (in list
    /// order). Returns the annotation that becomes the current one, `None`
    /// when none is left selected.
    pub fn click_annotation(&mut self, id: i64, ctrl: bool, shift: bool) -> Option<i64> {
        let index = |x: i64| self.annotations.iter().position(|a| a.id == x);
        if shift && let (Some(i), Some(j)) = (self.annotation_anchor.and_then(index), index(id)) {
            self.selected_annotations = self.annotations[i.min(j)..=i.max(j)].iter().map(|a| a.id).filter(|&x| x != id).collect();
        } else if ctrl && self.selected_annotations.contains(&id) {
            self.selected_annotations.retain(|&x| x != id);
            self.annotation_anchor = Some(id);
            return self.selected_annotation();
        } else {
            if !ctrl {
                self.selected_annotations.clear();
            }
            self.annotation_anchor = Some(id);
        }
        self.selected_annotations.retain(|&x| x != id);
        self.selected_annotations.push(id);
        Some(id)
    }

    /// Oligo of columns `[c0, c1)` in the given source, 5'→3' for the
    /// annotation type (reverse complement for reverse primers).
    pub fn oligo(&self, source: OligoSource, kind: crate::model::AnnotationKind, c0: usize, c1: usize) -> String {
        let t = self.template_row(source);
        let s = crate::seq::ungap(&t[c0.min(t.len())..c1.min(t.len())]);
        let s = if kind == crate::model::AnnotationKind::ReversePrimer { crate::seq::revcomp(&s) } else { s };
        String::from_utf8_lossy(&s).into_owned()
    }

    /// Row shown by a display item (group headers have none).
    pub fn item_row(&self, item: usize) -> Option<usize> {
        match self.items.get(item)? {
            Item::Seq { row, .. } => Some(*row),
            Item::Group { .. } => None,
        }
    }

    pub fn item_of_row(&self, row: usize) -> Option<usize> {
        self.items.iter().position(|it| matches!(it, Item::Seq { row: r, .. } if *r == row))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{AnnotationKind, DocKind, Row};

    #[test]
    fn annotation_clicks() {
        let doc = Document::new(vec![Row::from_gapped("s".into(), String::new(), b"ACGTACGTAC".to_vec())]);
        let info = DocInfo {
            id: 1,
            name: "d".into(),
            kind: DocKind::Sequences,
            n_rows: 1,
            width: doc.width,
            reference: None,
            created: String::new(),
            description: String::new(),
            folder: None,
        };
        // Ids 10..=14 in list order.
        let anns = (0..5)
            .map(|i| Annotation { id: 10 + i, name: format!("a{i}"), kind: AnnotationKind::Region, start: i as usize, end: i as usize + 2, note: String::new(), color: None })
            .collect();
        let mut od = OpenDoc::new(info, doc, anns, &ViewPrefs::default());
        // Shift without an anchor selects only the one clicked.
        assert_eq!(od.click_annotation(11, false, true), Some(11));
        assert_eq!(od.selected_annotations, [11]);
        // Ctrl adds; the one clicked becomes the current one.
        assert_eq!(od.click_annotation(13, true, false), Some(13));
        assert_eq!(od.selected_annotations, [11, 13]);
        // Shift selects from the anchor (13) to it, either way round.
        assert_eq!(od.click_annotation(14, false, true), Some(14));
        assert_eq!(od.selected_annotations, [13, 14]);
        assert_eq!(od.click_annotation(10, false, true), Some(10));
        assert_eq!(od.selected_annotations, [11, 12, 13, 10]);
        // Ctrl on a selected one removes it; the last one left is current.
        assert_eq!(od.click_annotation(10, true, false), Some(13));
        assert_eq!(od.selected_annotations, [11, 12, 13]);
        // A plain click selects it alone; Ctrl on it then leaves nothing.
        assert_eq!(od.click_annotation(12, false, false), Some(12));
        assert_eq!(od.selected_annotations, [12]);
        assert_eq!(od.click_annotation(12, true, false), None);
        assert!(od.selected_annotations.is_empty());
    }
}
