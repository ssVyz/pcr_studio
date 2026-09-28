//! State of the open document and the viewer.

use crate::display::{self, DisplayOptions, Item};
use crate::model::{Annotation, ColumnStats, ConsensusThreshold, DocInfo, Document, RefCoords};
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
        }
    }
}

pub const NAME_W_DEFAULT: f32 = 250.0;
pub const NAME_W_MIN: f32 = 120.0;
pub const NAME_W_MAX: f32 = 600.0;

/// Zoom level where letters are comfortable; shown as 100%.
pub const BASE_COL_W: f32 = 12.0;
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
    pub annotation_name: String,
    pub annotation_kind: crate::model::AnnotationKind,
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
            annotation_name: String::new(),
            annotation_kind: crate::model::AnnotationKind::ForwardPrimer,
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

    /// Top of display item `i` in content pixels (`i == items.len()` gives the total height).
    pub fn item_top(&self, i: usize) -> f32 {
        let g = self.group_prefix[i.min(self.items.len())] as f32;
        (i as f32 - g) * self.vp.row_h + g * self.group_h()
    }

    pub fn item_h(&self, i: usize) -> f32 {
        match self.items.get(i) {
            Some(Item::Group { .. }) => self.group_h(),
            _ => self.vp.row_h,
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
        self.vp.clamp(self.doc.width, h);
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
        let top = self.item_at_y(self.vp.scroll_y);
        self.vp.row_h = h;
        if let Some(i) = top {
            self.vp.scroll_y = self.item_top(i);
        }
        self.clamp_view();
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
        let rc = self.ref_coords.as_ref()?;
        let c1 = c1.min(rc.pos_at.len());
        if c1 == 0 || c0 >= c1 {
            return None;
        }
        let first = (c0..c1).find(|&c| rc.is_base[c]).map(|c| rc.pos_at[c])?;
        let last = rc.pos_at[(c1 - 1).min(rc.pos_at.len() - 1)];
        Some((first, last))
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
