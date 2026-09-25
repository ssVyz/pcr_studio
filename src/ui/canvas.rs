//! The alignment viewer: a canvas drawing only the visible part of the
//! alignment. Layout (top to bottom): overview strip, ruler, consensus,
//! identity graph, annotation lanes, pinned reference row, sequence rows.
//! A name column sits on the left, scrollbars on the right and bottom.

use super::Message;
use super::state::{Highlight, Hover, OpenDoc, ViewPrefs};
use super::style::Palette;
use crate::display::Item;
use crate::model::{Annotation, Row};
use crate::seq::GAP;
use iced::alignment::Vertical;
use iced::mouse::{self, Cursor, ScrollDelta};
use iced::widget::canvas::{self, Action, Event, Frame, Geometry, Path, Stroke, Text};
use iced::{Color, Font, Pixels, Point, Rectangle, Renderer, Size, Theme, keyboard};

pub const NAME_W: f32 = 250.0;
pub const SB: f32 = 12.0;
const OVERVIEW_H: f32 = 22.0;
const RULER_H: f32 = 24.0;
const IDENTITY_H: f32 = 42.0;
const LANE_H: f32 = 18.0;
const MAX_LANES: usize = 8;

#[derive(Debug, Clone)]
pub enum ViewMsg {
    Resized { canvas: (f32, f32), seq_w: f32, rows_h: f32 },
    Scroll { dx_cols: f32, dy_px: f32 },
    SetScroll { x: Option<f32>, y: Option<f32> },
    /// Zoom around `anchor_col`, which stays at `anchor_frac` of the sequence area.
    Zoom { factor: f32, anchor_col: f32, anchor_frac: f32 },
    Select { anchor: usize, to: usize, finished: bool },
    ClickItem { item: usize, ctrl: bool, shift: bool },
    ToggleGroup(String),
    Hover(Option<Hover>),
    ClickAnnotation(usize),
}

/// Vertical layout of the canvas for a given size.
#[derive(Debug, Clone)]
pub struct Lay {
    pub w: f32,
    pub h: f32,
    pub seq_x0: f32,
    pub seq_w: f32,
    pub overview: Option<(f32, f32)>,
    pub ruler: (f32, f32),
    pub consensus: Option<(f32, f32)>,
    pub identity: Option<(f32, f32)>,
    pub annotations: Option<(f32, f32)>,
    pub lanes: Vec<usize>,
    pub reference: Option<(f32, f32)>,
    pub rows_y0: f32,
    pub rows_h: f32,
}

/// Greedy lane assignment for overlapping annotations.
pub fn annotation_lanes(anns: &[Annotation]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..anns.len()).collect();
    order.sort_by_key(|&i| (anns[i].start, anns[i].end));
    let mut lane_end: Vec<usize> = Vec::new();
    let mut lanes = vec![0; anns.len()];
    for i in order {
        let a = &anns[i];
        let lane = match lane_end.iter().position(|&e| e <= a.start) {
            Some(l) => l,
            None if lane_end.len() < MAX_LANES => {
                lane_end.push(0);
                lane_end.len() - 1
            }
            None => MAX_LANES - 1,
        };
        lane_end[lane] = lane_end[lane].max(a.end);
        lanes[i] = lane;
    }
    lanes
}

pub fn layout(od: &OpenDoc, prefs: &ViewPrefs, size: Size) -> Lay {
    let seq_x0 = NAME_W;
    let seq_w = (size.width - NAME_W - SB).max(10.0);
    let mut y = 0.0;
    let mut take = |h: f32| {
        let r = (y, h);
        y += h;
        r
    };
    let overview = prefs.show_overview.then(|| take(OVERVIEW_H));
    let ruler = take(RULER_H);
    let track_h = od.vp.row_h.clamp(14.0, 20.0);
    let consensus = prefs.show_consensus.then(|| take(track_h));
    let identity = prefs.show_identity.then(|| take(IDENTITY_H));
    let lanes = annotation_lanes(&od.annotations);
    let n_lanes = lanes.iter().map(|l| l + 1).max().unwrap_or(0);
    let annotations = (prefs.show_annotations && n_lanes > 0).then(|| take(n_lanes as f32 * LANE_H + 4.0));
    let reference = od.reference.map(|_| take(track_h.max(od.vp.row_h)));
    let rows_y0 = y + 2.0;
    let rows_h = (size.height - SB - rows_y0).max(10.0);
    Lay { w: size.width, h: size.height, seq_x0, seq_w, overview, ruler, consensus, identity, annotations, lanes, reference, rows_y0, rows_h }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
enum Drag {
    #[default]
    None,
    Select {
        anchor: usize,
    },
    VThumb {
        grab: f32,
    },
    HThumb {
        grab: f32,
    },
    Overview,
    Pan {
        x: f32,
        y: f32,
        sx: f32,
        sy: f32,
    },
}

#[derive(Default)]
pub struct CanvasState {
    drag: Drag,
    modifiers: keyboard::Modifiers,
    last_hover: Option<Hover>,
}

pub struct AlignmentView<'a> {
    pub od: &'a OpenDoc,
    pub prefs: &'a ViewPrefs,
    pub pal: Palette,
}

impl AlignmentView<'_> {
    fn col_at(&self, lay: &Lay, x: f32) -> f32 {
        self.od.vp.scroll_x + (x - lay.seq_x0) / self.od.vp.col_w
    }

    fn col_index(&self, lay: &Lay, x: f32) -> usize {
        let w = self.od.width();
        (self.col_at(lay, x).floor().max(0.0) as usize).min(w.saturating_sub(1))
    }

    fn x_of(&self, lay: &Lay, col: f32) -> f32 {
        lay.seq_x0 + (col - self.od.vp.scroll_x) * self.od.vp.col_w
    }

    fn item_at(&self, lay: &Lay, y: f32) -> Option<usize> {
        if y < lay.rows_y0 || y > lay.rows_y0 + lay.rows_h {
            return None;
        }
        self.od.item_at_y(self.od.vp.scroll_y + (y - lay.rows_y0))
    }

    fn vthumb(&self, lay: &Lay) -> (f32, f32) {
        let total = self.od.content_h().max(1.0);
        let vis = lay.rows_h.min(total);
        let th = (lay.rows_h * vis / total).max(24.0).min(lay.rows_h);
        let max_scroll = (total - lay.rows_h + self.od.vp.row_h * 0.5).max(0.0);
        let frac = if max_scroll > 0.0 { self.od.vp.scroll_y / max_scroll } else { 0.0 };
        (lay.rows_y0 + (lay.rows_h - th) * frac.clamp(0.0, 1.0), th)
    }

    fn hthumb(&self, lay: &Lay) -> (f32, f32) {
        let total = self.od.width().max(1) as f32;
        let vis = (lay.seq_w / self.od.vp.col_w).min(total);
        let tw = (lay.seq_w * vis / total).max(24.0).min(lay.seq_w);
        let max_scroll = (total - vis).max(0.0);
        let frac = if max_scroll > 0.0 { self.od.vp.scroll_x / max_scroll } else { 0.0 };
        (lay.seq_x0 + (lay.seq_w - tw) * frac.clamp(0.0, 1.0), tw)
    }

    fn hover_at(&self, lay: &Lay, p: Point) -> Option<Hover> {
        if p.x < lay.seq_x0 || p.x > lay.seq_x0 + lay.seq_w || self.od.width() == 0 {
            return None;
        }
        let col = self.col_index(lay, p.x);
        let on_reference = lay.reference.is_some_and(|(y, h)| p.y >= y && p.y < y + h);
        Some(Hover { col, item: self.item_at(lay, p.y), on_reference })
    }
}

fn in_band(y: f32, band: Option<(f32, f32)>) -> bool {
    band.is_some_and(|(y0, h)| y >= y0 && y < y0 + h)
}

impl canvas::Program<Message> for AlignmentView<'_> {
    type State = CanvasState;

    fn update(&self, state: &mut CanvasState, event: &Event, bounds: Rectangle, cursor: Cursor) -> Option<Action<Message>> {
        let size = bounds.size();
        let lay = layout(self.od, self.prefs, size);
        if (size.width, size.height) != self.od.vp.canvas || (lay.seq_w - self.od.vp.seq_w).abs() > 0.5 || (lay.rows_h - self.od.vp.rows_h).abs() > 0.5 {
            return Some(Action::publish(Message::View(ViewMsg::Resized {
                canvas: (size.width, size.height),
                seq_w: lay.seq_w,
                rows_h: lay.rows_h,
            })));
        }
        if let Event::Keyboard(keyboard::Event::ModifiersChanged(m)) = event {
            state.modifiers = *m;
            return None;
        }
        let pos = cursor.position_in(bounds);
        let vp = &self.od.vp;
        let msg = |m: ViewMsg| Some(Action::publish(Message::View(m)).and_capture());
        match event {
            Event::Mouse(mouse::Event::WheelScrolled { delta }) => {
                let p = pos?;
                let (dx, dy) = match *delta {
                    ScrollDelta::Lines { x, y } => (x * 40.0, y * 40.0),
                    ScrollDelta::Pixels { x, y } => (x, y),
                };
                if state.modifiers.control() {
                    let factor = 1.2f32.powf(dy / 40.0);
                    let x = p.x.clamp(lay.seq_x0, lay.seq_x0 + lay.seq_w);
                    return msg(ViewMsg::Zoom {
                        factor,
                        anchor_col: self.col_at(&lay, x),
                        anchor_frac: (x - lay.seq_x0) / lay.seq_w,
                    });
                }
                let (dx, dy) = if state.modifiers.shift() && dx == 0.0 { (dy, 0.0) } else { (dx, dy) };
                return msg(ViewMsg::Scroll { dx_cols: -dx * 2.0 / vp.col_w, dy_px: -dy * 1.5 });
            }
            Event::Mouse(mouse::Event::ButtonPressed(mouse::Button::Middle)) => {
                let p = pos?;
                state.drag = Drag::Pan { x: p.x, y: p.y, sx: vp.scroll_x, sy: vp.scroll_y };
                return Some(Action::capture());
            }
            Event::Mouse(mouse::Event::ButtonPressed(mouse::Button::Left)) => {
                let p = pos?;
                let width = self.od.width();
                // Vertical scrollbar.
                if p.x >= lay.w - SB && p.y >= lay.rows_y0 {
                    let (ty, th) = self.vthumb(&lay);
                    if p.y >= ty && p.y <= ty + th {
                        state.drag = Drag::VThumb { grab: p.y - ty };
                        return Some(Action::capture());
                    }
                    let page = lay.rows_h * 0.9;
                    return msg(ViewMsg::Scroll { dx_cols: 0.0, dy_px: if p.y < ty { -page } else { page } });
                }
                // Horizontal scrollbar.
                if p.y >= lay.h - SB && p.x >= lay.seq_x0 {
                    let (tx, tw) = self.hthumb(&lay);
                    if p.x >= tx && p.x <= tx + tw {
                        state.drag = Drag::HThumb { grab: p.x - tx };
                        return Some(Action::capture());
                    }
                    let page = lay.seq_w / vp.col_w * 0.9;
                    return msg(ViewMsg::Scroll { dx_cols: if p.x < tx { -page } else { page }, dy_px: 0.0 });
                }
                if in_band(p.y, lay.overview) && p.x >= lay.seq_x0 {
                    state.drag = Drag::Overview;
                    let col = (p.x - lay.seq_x0) / lay.seq_w * width as f32;
                    return msg(ViewMsg::SetScroll { x: Some(col - lay.seq_w / vp.col_w / 2.0), y: None });
                }
                if p.x < lay.seq_x0 {
                    if let Some(item) = self.item_at(&lay, p.y) {
                        if let Some(Item::Group { label, .. }) = self.od.items.get(item) {
                            return msg(ViewMsg::ToggleGroup(label.clone()));
                        }
                        let m = state.modifiers;
                        return msg(ViewMsg::ClickItem { item, ctrl: m.control(), shift: m.shift() });
                    }
                    return None;
                }
                if p.x > lay.seq_x0 + lay.seq_w || width == 0 {
                    return None;
                }
                if let Some((y0, _)) = lay.annotations.filter(|_| in_band(p.y, lay.annotations)) {
                    let lane = ((p.y - y0 - 2.0) / LANE_H).floor().max(0.0) as usize;
                    let col = self.col_index(&lay, p.x);
                    if let Some(i) = self.od.annotations.iter().enumerate().position(|(i, a)| lay.lanes[i] == lane && col >= a.start && col < a.end) {
                        return msg(ViewMsg::ClickAnnotation(i));
                    }
                }
                let col = self.col_index(&lay, p.x);
                let anchor = match (state.modifiers.shift(), self.od.selection) {
                    (true, Some(sel)) => sel.anchor,
                    _ => col,
                };
                state.drag = Drag::Select { anchor };
                return msg(ViewMsg::Select { anchor, to: col, finished: false });
            }
            Event::Mouse(mouse::Event::CursorMoved { .. }) => {
                let Some(p) = cursor.position_in(bounds).or_else(|| cursor.position().map(|q| Point::new(q.x - bounds.x, q.y - bounds.y)))
                else {
                    return None;
                };
                match state.drag {
                    Drag::Select { anchor } => {
                        let col = self.col_index(&lay, p.x.clamp(lay.seq_x0, lay.seq_x0 + lay.seq_w - 1.0));
                        return msg(ViewMsg::Select { anchor, to: col, finished: false });
                    }
                    Drag::VThumb { grab } => {
                        let (_, th) = self.vthumb(&lay);
                        let max_scroll = (self.od.content_h() - lay.rows_h + vp.row_h * 0.5).max(0.0);
                        let span = (lay.rows_h - th).max(1.0);
                        let frac = ((p.y - grab - lay.rows_y0) / span).clamp(0.0, 1.0);
                        return msg(ViewMsg::SetScroll { x: None, y: Some(frac * max_scroll) });
                    }
                    Drag::HThumb { grab } => {
                        let (_, tw) = self.hthumb(&lay);
                        let total = self.od.width() as f32;
                        let vis = lay.seq_w / vp.col_w;
                        let span = (lay.seq_w - tw).max(1.0);
                        let frac = ((p.x - grab - lay.seq_x0) / span).clamp(0.0, 1.0);
                        return msg(ViewMsg::SetScroll { x: Some(frac * (total - vis).max(0.0)), y: None });
                    }
                    Drag::Overview => {
                        let col = (p.x - lay.seq_x0) / lay.seq_w * self.od.width() as f32;
                        return msg(ViewMsg::SetScroll { x: Some(col - lay.seq_w / vp.col_w / 2.0), y: None });
                    }
                    Drag::Pan { x, y, sx, sy } => {
                        return msg(ViewMsg::SetScroll {
                            x: Some(sx - (p.x - x) / vp.col_w),
                            y: Some(sy - (p.y - y)),
                        });
                    }
                    Drag::None => {
                        let hover = if cursor.is_over(bounds) { self.hover_at(&lay, p) } else { None };
                        if hover != state.last_hover {
                            state.last_hover = hover;
                            return Some(Action::publish(Message::View(ViewMsg::Hover(hover))));
                        }
                    }
                }
            }
            Event::Mouse(mouse::Event::ButtonReleased(_)) => {
                let drag = std::mem::take(&mut state.drag);
                if let Drag::Select { anchor } = drag {
                    let p = cursor.position().map(|q| Point::new(q.x - bounds.x, q.y - bounds.y)).unwrap_or(Point::ORIGIN);
                    let col = self.col_index(&lay, p.x.clamp(lay.seq_x0, lay.seq_x0 + lay.seq_w - 1.0));
                    return msg(ViewMsg::Select { anchor, to: col, finished: true });
                }
                if drag != Drag::None {
                    return Some(Action::capture());
                }
            }
            Event::Mouse(mouse::Event::CursorLeft)
                if state.last_hover.is_some() => {
                    state.last_hover = None;
                    return Some(Action::publish(Message::View(ViewMsg::Hover(None))));
                }
            _ => {}
        }
        None
    }

    fn draw(&self, _state: &CanvasState, renderer: &Renderer, _theme: &Theme, bounds: Rectangle, _cursor: Cursor) -> Vec<Geometry> {
        let mut f = Frame::new(renderer, bounds.size());
        let lay = layout(self.od, self.prefs, bounds.size());
        let pal = &self.pal;
        let full = Rectangle::with_size(bounds.size());
        // Everything goes through clipped sub-frames: canvas text renders above
        // all shapes, so letters must be clipped to their region, and pasted
        // sub-frames keep their drawing order (shapes drawn directly on `f`
        // would land on top of them).
        f.with_clip(full, |f| f.fill_rectangle(Point::ORIGIN, bounds.size(), pal.bg));
        if self.od.width() == 0 {
            return vec![f.into_geometry()];
        }
        let d = Drawer { v: self, lay: &lay, pal };
        let seq_clip = |y: f32, h: f32| Rectangle { x: lay.seq_x0, y, width: lay.seq_w, height: h.max(0.0) };
        f.with_clip(seq_clip(lay.rows_y0, lay.rows_h), |f| d.rows(f));
        f.with_clip(seq_clip(0.0, lay.rows_y0), |f| d.headers(f));
        f.with_clip(seq_clip(0.0, lay.h - SB), |f| d.selection(f));
        f.with_clip(Rectangle { x: 0.0, y: 0.0, width: lay.seq_x0, height: lay.h - SB }, |f| d.names(f));
        f.with_clip(full, |f| d.scrollbars(f));
        vec![f.into_geometry()]
    }

    fn mouse_interaction(&self, state: &CanvasState, bounds: Rectangle, cursor: Cursor) -> mouse::Interaction {
        match state.drag {
            Drag::Pan { .. } => return mouse::Interaction::Grabbing,
            Drag::Select { .. } => return mouse::Interaction::Text,
            _ => {}
        }
        let Some(p) = cursor.position_in(bounds) else { return mouse::Interaction::default() };
        let lay = layout(self.od, self.prefs, bounds.size());
        if p.x < lay.seq_x0 {
            if self.item_at(&lay, p.y).is_some() {
                return mouse::Interaction::Pointer;
            }
        } else if p.x < lay.seq_x0 + lay.seq_w && p.y < lay.h - SB {
            if in_band(p.y, lay.overview) || in_band(p.y, lay.annotations) {
                return mouse::Interaction::Pointer;
            }
            return mouse::Interaction::Crosshair;
        }
        mouse::Interaction::default()
    }
}

/// How one alignment cell is drawn.
#[derive(Clone, Copy, PartialEq)]
enum Cell {
    Blank,
    /// Colored block with a letter.
    Base(u8),
    /// Same as the comparison sequence.
    Same(u8),
    Gap,
    SameGap,
}

struct Drawer<'a, 'b> {
    v: &'a AlignmentView<'b>,
    lay: &'a Lay,
    pal: &'a Palette,
}

/// Merges horizontally adjacent rectangles of one color.
struct Runs {
    cur: Option<(f32, f32, Color)>,
    y: f32,
    h: f32,
}

impl Runs {
    fn new(y: f32, h: f32) -> Runs {
        Runs { cur: None, y, h }
    }
    fn push(&mut self, f: &mut Frame, x0: f32, x1: f32, c: Color) {
        match &mut self.cur {
            Some((_, e, cc)) if *cc == c && (x0 - *e).abs() < 0.01 => *e = x1,
            _ => {
                self.flush(f);
                self.cur = Some((x0, x1, c));
            }
        }
    }
    fn flush(&mut self, f: &mut Frame) {
        if let Some((a, b, c)) = self.cur.take() {
            f.fill_rectangle(Point::new(a, self.y), Size::new((b - a).max(0.5), self.h), c);
        }
    }
}

impl Drawer<'_, '_> {
    fn od(&self) -> &OpenDoc {
        self.v.od
    }

    /// Visible display items as (index, screen y, height).
    fn visible_items(&self) -> Vec<(usize, f32, f32)> {
        let od = self.od();
        let lay = self.lay;
        let mut out = Vec::new();
        let Some(mut i) = od.item_at_y(od.vp.scroll_y.max(0.0)) else { return out };
        while i < od.items.len() {
            let y = lay.rows_y0 + od.item_top(i) - od.vp.scroll_y;
            if y > lay.rows_y0 + lay.rows_h {
                break;
            }
            out.push((i, y, od.item_h(i)));
            i += 1;
        }
        out
    }

    fn visible_cols(&self) -> (usize, usize) {
        let vp = &self.od().vp;
        let c0 = vp.scroll_x.floor().max(0.0) as usize;
        let c1 = ((vp.scroll_x + self.lay.seq_w / vp.col_w).ceil() as usize + 1).min(self.od().width());
        (c0.min(c1), c1)
    }

    fn classify(&self, b: Option<u8>, cmp: Option<u8>, highlight: bool) -> Cell {
        match b {
            None => Cell::Blank,
            Some(GAP) => {
                if highlight && cmp == Some(GAP) {
                    Cell::SameGap
                } else {
                    Cell::Gap
                }
            }
            Some(x) => {
                if highlight && cmp == Some(x) {
                    Cell::Same(x)
                } else {
                    Cell::Base(x)
                }
            }
        }
    }

    /// Background of a cell; `letters` tells whether glyphs are drawn on top.
    fn cell_fill(&self, cell: Cell, letters: bool) -> Option<Color> {
        let prefs = self.v.prefs;
        let pal = self.pal;
        match cell {
            Cell::Blank => None,
            Cell::Base(b) => Some(if prefs.color_bases { pal.base(b) } else { pal.same }),
            // Dots replace the fill only when they are visible; otherwise a faint
            // fill keeps covered regions distinguishable from missing coverage.
            Cell::Same(_) => (!prefs.use_dots || !letters).then_some(pal.same),
            Cell::Gap => prefs.highlight_gaps.then_some(pal.gap_strong),
            Cell::SameGap => None,
        }
    }

    fn letter(&self, f: &mut Frame, x: f32, y: f32, h: f32, ch: char, color: Color, size: f32) {
        f.fill_text(Text {
            content: ch.to_string(),
            position: Point::new(x, y + h / 2.0),
            color,
            size: Pixels(size),
            font: Font::MONOSPACE,
            align_x: iced::widget::text::Alignment::Center,
            align_y: Vertical::Center,
            ..Text::default()
        });
    }

    fn label(&self, f: &mut Frame, x: f32, y: f32, h: f32, s: &str, color: Color, size: f32, bold: bool) {
        let font = if bold { Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT } } else { Font::DEFAULT };
        f.fill_text(Text {
            content: s.to_string(),
            position: Point::new(x, y + h / 2.0),
            color,
            size: Pixels(size),
            font,
            align_y: Vertical::Center,
            ..Text::default()
        });
    }

    /// Draws one sequence track (a row, the reference or the consensus).
    fn track(&self, f: &mut Frame, y: f32, h: f32, at: &dyn Fn(usize) -> Option<u8>, cmp: Option<&dyn Fn(usize) -> Option<u8>>) {
        let vp = &self.od().vp;
        let col_w = vp.col_w;
        let (c0, c1) = self.visible_cols();
        let highlight = cmp.is_some();
        let pad = if h >= 8.0 { 1.0 } else { 0.0 };
        let (by, bh) = (y + pad / 2.0, h - pad);
        let mut runs = Runs::new(by, bh);
        let gap_line_h = (h * 0.12).clamp(1.0, 2.0);
        let mut gaps = Runs::new(y + h / 2.0 - gap_line_h / 2.0, gap_line_h);
        if col_w >= 1.0 {
            let letters = col_w >= 7.0 && h >= 9.0;
            let size = (col_w.min(h) * 0.82).clamp(7.0, 15.0);
            let mut glyphs: Vec<(f32, char, Color)> = Vec::new();
            for c in c0..c1 {
                let b = at(c);
                let cell = self.classify(b, cmp.and_then(|g| g(c)), highlight);
                let x = self.v.x_of(self.lay, c as f32);
                if let Some(fill) = self.cell_fill(cell, letters) {
                    runs.push(f, x, x + col_w, fill);
                } else {
                    runs.flush(f);
                }
                match cell {
                    Cell::Gap if !self.v.prefs.highlight_gaps => gaps.push(f, x, x + col_w, self.pal.gap),
                    Cell::SameGap if !self.v.prefs.use_dots => gaps.push(f, x, x + col_w, self.pal.gap),
                    _ => gaps.flush(f),
                }
                if letters {
                    match cell {
                        Cell::Base(b) => glyphs.push((x + col_w / 2.0, b as char, self.pal.base_text)),
                        Cell::Same(b) => {
                            let ch = if self.v.prefs.use_dots { '.' } else { b as char };
                            glyphs.push((x + col_w / 2.0, ch, self.pal.same_text));
                        }
                        Cell::Gap if self.v.prefs.highlight_gaps => glyphs.push((x + col_w / 2.0, '-', self.pal.bg)),
                        _ => {}
                    }
                }
            }
            runs.flush(f);
            gaps.flush(f);
            for (x, ch, color) in glyphs {
                self.letter(f, x, y, h, ch, color, size);
            }
        } else {
            // Several columns per pixel: summarize each pixel column.
            let px_n = self.lay.seq_w.ceil() as usize;
            let width = self.od().width();
            for px in 0..px_n {
                let x = self.lay.seq_x0 + px as f32;
                let ca = self.v.col_at(self.lay, x).floor().max(0.0) as usize;
                let cb = (self.v.col_at(self.lay, x + 1.0).floor().max(0.0) as usize).max(ca + 1).min(width);
                if ca >= width {
                    break;
                }
                let mut color: Option<Color> = None;
                let mut covered = false;
                let mut gap_seen = false;
                if highlight {
                    let g = cmp.unwrap();
                    for c in ca..cb {
                        let Some(b) = at(c) else { continue };
                        covered = true;
                        let r = g(c);
                        if r != Some(b) {
                            if b == GAP {
                                gap_seen = true;
                            } else {
                                color = Some(if self.v.prefs.color_bases { self.pal.base(b) } else { self.pal.text_muted });
                                break;
                            }
                        }
                    }
                } else {
                    let mid = (ca + cb) / 2;
                    match at(mid) {
                        Some(GAP) => {
                            covered = true;
                            gap_seen = true;
                        }
                        Some(b) => {
                            covered = true;
                            color = Some(if self.v.prefs.color_bases { self.pal.base(b) } else { self.pal.same });
                        }
                        None => {
                            covered = (ca..cb).any(|c| at(c).is_some());
                        }
                    }
                }
                let color = match (color, gap_seen, covered) {
                    (Some(c), _, _) => Some(c),
                    (None, true, _) => Some(if self.v.prefs.highlight_gaps { self.pal.gap_strong } else { self.pal.gap }),
                    (None, false, true) => Some(self.pal.same),
                    _ => None,
                };
                match color {
                    Some(c) => runs.push(f, x, x + 1.0, c),
                    None => runs.flush(f),
                }
            }
            runs.flush(f);
        }
    }

    fn rows(&self, f: &mut Frame) {
        let od = self.od();
        let vp = &od.vp;
        let lay = self.lay;
        let doc = &od.doc;
        let reference = od.reference.map(|r| &doc.rows[r]);
        let consensus = &od.consensus;
        let cons_at = |c: usize| -> Option<u8> { consensus.get(c).copied().filter(|&b| b != b' ') };
        let ref_at = |c: usize| -> Option<u8> { reference.and_then(|r| r.at(c)) };
        let cmp: Option<&dyn Fn(usize) -> Option<u8>> = match (self.v.prefs.highlight, reference) {
            (Highlight::Reference, Some(_)) => Some(&ref_at),
            (Highlight::Consensus, _) => Some(&cons_at),
            _ => None,
        };
        for (i, y, h) in self.visible_items() {
            match &od.items[i] {
                Item::Group { .. } => {
                    f.fill_rectangle(Point::new(lay.seq_x0, y), Size::new(lay.seq_w, h), self.pal.group_bg);
                }
                Item::Seq { row, .. } => {
                    let r: &Row = &doc.rows[*row];
                    if od.selected_rows.contains(row) {
                        f.fill_rectangle(Point::new(lay.seq_x0, y), Size::new(lay.seq_w, vp.row_h), self.pal.row_selected);
                    }
                    self.track(f, y, vp.row_h, &|c| r.at(c), cmp);
                }
            }
        }
        // Current search hit.
        if let Some(h) = od.search.current.and_then(|i| od.search.hits.get(i))
            && let Some(item) = h.row.and_then(|r| od.item_of_row(r)) {
                let y = lay.rows_y0 + od.item_top(item) - vp.scroll_y;
                if y >= lay.rows_y0 - vp.row_h && y <= lay.rows_y0 + lay.rows_h {
                    self.hit_box(f, h.c0, h.c1, y, vp.row_h);
                }
            }
    }

    fn hit_box(&self, f: &mut Frame, c0: usize, c1: usize, y: f32, h: f32) {
        let x0 = self.v.x_of(self.lay, c0 as f32);
        let x1 = self.v.x_of(self.lay, c1 as f32).max(x0 + 2.0);
        f.stroke(&Path::rectangle(Point::new(x0, y), Size::new(x1 - x0, h)), Stroke::default().with_color(self.pal.hit).with_width(2.0));
    }

    fn headers(&self, f: &mut Frame) {
        let od = self.od();
        let lay = self.lay;
        let pal = self.pal;
        let vp = &od.vp;
        let width = od.width();
        let header_h = lay.rows_y0;
        f.fill_rectangle(Point::new(lay.seq_x0, 0.0), Size::new(lay.seq_w + SB, header_h), pal.header_bg);
        let (c0, c1) = self.visible_cols();

        // Overview strip: identity profile over the whole alignment.
        if let Some((y, h)) = lay.overview {
            let bar_y = y + 4.0;
            let bar_h = h - 8.0;
            let px_n = lay.seq_w.floor() as usize;
            let mode = self.v.prefs.graph;
            let track = od.graph_track(mode);
            let mut runs = Runs::new(bar_y, bar_h);
            for px in 0..px_n {
                let ca = px * width / px_n.max(1);
                let cb = ((px + 1) * width / px_n.max(1)).max(ca + 1).min(width);
                let v = bin_value(&track[ca..cb], mode);
                let x = lay.seq_x0 + px as f32;
                match v {
                    Some(v) => runs.push(f, x, x + 1.0, pal.graph(mode, v)),
                    None => runs.push(f, x, x + 1.0, pal.grid),
                }
            }
            runs.flush(f);
            let to_x = |c: f32| lay.seq_x0 + c / width as f32 * lay.seq_w;
            for a in &od.annotations {
                let x0 = to_x(a.start as f32);
                let x1 = to_x(a.end as f32).max(x0 + 2.0);
                f.fill_rectangle(Point::new(x0, y), Size::new(x1 - x0, 3.0), pal.annotation(a.kind));
            }
            for h in od.search.hits.iter().take(2000) {
                f.fill_rectangle(Point::new(to_x(h.c0 as f32), y + h_off(h.row.is_none())), Size::new(1.5, 3.0), pal.hit);
            }
            if let Some(sel) = od.selection {
                let x0 = to_x(sel.c0 as f32);
                f.fill_rectangle(Point::new(x0 - 1.0, y + 1.0), Size::new((to_x(sel.c1 as f32) - x0).max(2.0) + 1.0, h - 2.0), pal.selection);
            }
            let vx0 = to_x(vp.scroll_x);
            let vx1 = to_x((vp.scroll_x + lay.seq_w / vp.col_w).min(width as f32));
            f.stroke(
                &Path::rectangle(Point::new(vx0, y + 2.0), Size::new((vx1 - vx0).max(3.0), h - 4.0)),
                Stroke::default().with_color(pal.text).with_width(1.5),
            );
        }

        // Ruler.
        let (ry, rh) = lay.ruler;
        f.fill_rectangle(Point::new(lay.seq_x0, ry + rh - 1.0), Size::new(lay.seq_w, 1.0), pal.grid);
        let step = nice_step(80.0 / vp.col_w);
        let ref_coords = od.ref_coords.as_ref();
        let sel_xs: Vec<f32> = od
            .selection
            .map(|s| vec![self.v.x_of(lay, s.c0 as f32 + 0.5), self.v.x_of(lay, (s.c1 - 1) as f32 + 0.5)])
            .unwrap_or_default();
        let mut last_label_x = f32::NEG_INFINITY;
        for c in c0..c1 {
            let (num, show) = match ref_coords {
                Some(rc) => (rc.pos_at[c] as usize, rc.is_base[c] && ((rc.pos_at[c] as usize).is_multiple_of(step) || rc.pos_at[c] == 1)),
                None => (c + 1, (c + 1) % step == 0 || c == 0),
            };
            if !show {
                continue;
            }
            let x = self.v.x_of(lay, c as f32 + 0.5);
            f.fill_rectangle(Point::new(x, ry + rh - 6.0), Size::new(1.0, 5.0), pal.text_muted);
            if x - last_label_x > 60.0 && sel_xs.iter().all(|sx| (sx - x).abs() > 44.0) {
                f.fill_text(Text {
                    content: group_digits(num),
                    position: Point::new(x, ry + 3.0),
                    color: pal.text_muted,
                    size: Pixels(11.0),
                    align_x: iced::widget::text::Alignment::Center,
                    ..Text::default()
                });
                last_label_x = x;
            }
        }
        // Selection bounds on the ruler.
        if let Some(sel) = od.selection {
            for c in [sel.c0, sel.c1 - 1] {
                let num = match ref_coords {
                    Some(rc) => rc.pos_at[c] as usize,
                    None => c + 1,
                };
                let x = self.v.x_of(lay, c as f32 + 0.5);
                if x >= lay.seq_x0 && x <= lay.seq_x0 + lay.seq_w {
                    f.fill_text(Text {
                        content: group_digits(num),
                        position: Point::new(x, ry + 2.0),
                        color: pal.selection_edge,
                        size: Pixels(11.0),
                        font: Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT },
                        align_x: iced::widget::text::Alignment::Center,
                        ..Text::default()
                    });
                }
            }
        }

        // Consensus.
        if let Some((y, h)) = lay.consensus {
            f.fill_rectangle(Point::new(lay.seq_x0, y), Size::new(lay.seq_w, h), pal.consensus_bg);
            let cons = &od.consensus;
            self.track(f, y + 1.0, h - 2.0, &|c| cons.get(c).copied().filter(|&b| b != b' '), None);
            if let Some(h_) = od.search.current.and_then(|i| od.search.hits.get(i)).filter(|h| h.row.is_none()) {
                self.hit_box(f, h_.c0, h_.c1, y, h);
            }
        }

        // Identity / conservation graph.
        if let Some((y, h)) = lay.identity {
            let base = y + h - 2.0;
            let full = h - 6.0;
            let mode = self.v.prefs.graph;
            let track = od.graph_track(mode);
            let bar = |f: &mut Frame, x: f32, w: f32, v: f32| {
                let bh = (full * graph_height(mode, v)).max(1.0);
                f.fill_rectangle(Point::new(x, base - bh), Size::new(w, bh), pal.graph(mode, v));
            };
            if vp.col_w >= 1.0 {
                for c in c0..c1 {
                    if let Some(v) = track[c] {
                        bar(f, self.v.x_of(lay, c as f32), vp.col_w, v);
                    }
                }
            } else {
                let px_n = lay.seq_w.ceil() as usize;
                for px in 0..px_n {
                    let x = lay.seq_x0 + px as f32;
                    let ca = self.v.col_at(lay, x).floor().max(0.0) as usize;
                    let cb = (self.v.col_at(lay, x + 1.0).floor().max(0.0) as usize).max(ca + 1).min(width);
                    if ca >= width {
                        break;
                    }
                    // Zoomed out: the weakest column of the pixel, so variable sites stay visible.
                    if let Some(v) = track[ca..cb].iter().flatten().copied().reduce(f32::min) {
                        bar(f, x, 1.0, v);
                    }
                }
            }
            f.fill_rectangle(Point::new(lay.seq_x0, y + h - 1.0), Size::new(lay.seq_w, 1.0), pal.grid);
        }

        // Annotations.
        if let Some((y0, _)) = lay.annotations {
            for (i, a) in od.annotations.iter().enumerate() {
                if a.end <= c0 || a.start >= c1 {
                    continue;
                }
                let y = y0 + 2.0 + lay.lanes[i] as f32 * LANE_H;
                let x0 = self.v.x_of(lay, a.start as f32).max(lay.seq_x0 - 20.0);
                let x1 = self.v.x_of(lay, a.end as f32).min(lay.seq_x0 + lay.seq_w + 20.0);
                self.arrow(f, x0, x1, y + 1.0, LANE_H - 3.0, a);
            }
        }

        // Pinned reference row.
        if let (Some((y, h)), Some(r)) = (lay.reference, od.reference) {
            f.fill_rectangle(Point::new(lay.seq_x0, y), Size::new(lay.seq_w, h), pal.reference_bg);
            let row = &od.doc.rows[r];
            self.track(f, y + 1.0, h - 2.0, &|c| row.at(c), None);
            if let Some(h_) = od.search.current.and_then(|i| od.search.hits.get(i)).filter(|h| h.row == Some(r)) {
                self.hit_box(f, h_.c0, h_.c1, y, h);
            }
            f.fill_rectangle(Point::new(lay.seq_x0, y + h), Size::new(lay.seq_w, 1.0), pal.grid);
        }
    }

    fn arrow(&self, f: &mut Frame, x0: f32, x1: f32, y: f32, h: f32, a: &Annotation) {
        let color = self.pal.annotation(a.kind);
        let head = (h * 0.6).min((x1 - x0) * 0.4).max(0.0);
        let path = Path::new(|p| match a.kind.direction() {
            1 => {
                p.move_to(Point::new(x0, y));
                p.line_to(Point::new(x1 - head, y));
                p.line_to(Point::new(x1, y + h / 2.0));
                p.line_to(Point::new(x1 - head, y + h));
                p.line_to(Point::new(x0, y + h));
                p.close();
            }
            -1 => {
                p.move_to(Point::new(x1, y));
                p.line_to(Point::new(x0 + head, y));
                p.line_to(Point::new(x0, y + h / 2.0));
                p.line_to(Point::new(x0 + head, y + h));
                p.line_to(Point::new(x1, y + h));
                p.close();
            }
            _ => p.rectangle(Point::new(x0, y), Size::new(x1 - x0, h)),
        });
        f.fill(&path, color);
        let text_w = a.name.chars().count() as f32 * 6.5;
        if x1 - x0 > text_w + 8.0 {
            f.fill_text(Text {
                content: a.name.clone(),
                position: Point::new((x0 + x1) / 2.0, y + h / 2.0),
                color: Color::WHITE,
                size: Pixels(11.0),
                align_x: iced::widget::text::Alignment::Center,
                align_y: Vertical::Center,
                ..Text::default()
            });
        }
    }

    fn selection(&self, f: &mut Frame) {
        let od = self.od();
        let lay = self.lay;
        let Some(sel) = od.selection else { return };
        let x0 = self.v.x_of(lay, sel.c0 as f32).max(lay.seq_x0 - 2.0);
        let x1 = self.v.x_of(lay, sel.c1 as f32).min(lay.seq_x0 + lay.seq_w + 2.0);
        if x1 < lay.seq_x0 || x0 > lay.seq_x0 + lay.seq_w {
            return;
        }
        let y0 = lay.ruler.0 + lay.ruler.1;
        let y1 = lay.h - SB;
        let w = (x1 - x0).max(1.5);
        f.fill_rectangle(Point::new(x0, y0), Size::new(w, y1 - y0), self.pal.selection);
        f.fill_rectangle(Point::new(x0, y0), Size::new(1.0, y1 - y0), self.pal.selection_edge);
        f.fill_rectangle(Point::new(x0 + w - 1.0, y0), Size::new(1.0, y1 - y0), self.pal.selection_edge);
    }

    fn names(&self, f: &mut Frame) {
        let od = self.od();
        let lay = self.lay;
        let pal = self.pal;
        let vp = &od.vp;
        f.fill_rectangle(Point::ORIGIN, Size::new(lay.seq_x0, lay.h), pal.names_bg);
        f.fill_rectangle(Point::new(lay.seq_x0 - 1.0, 0.0), Size::new(1.0, lay.h), pal.grid);
        let x = 10.0;
        let head = |f: &mut Frame, band: Option<(f32, f32)>, s: &str| {
            if let Some((y, h)) = band {
                self.label(f, x, y, h, s, pal.text_muted, 12.0, false);
            }
        };
        head(f, lay.overview, "Overview");
        head(f, Some(lay.ruler), if od.reference.is_some() { "Reference position" } else { "Column" });
        head(f, lay.consensus, "Consensus");
        head(f, lay.identity, match self.v.prefs.graph {
            super::state::GraphMode::Conservation => "Conservation",
            super::state::GraphMode::Identity => "Identity",
        });
        head(f, lay.annotations, "Annotations");
        if let (Some((y, h)), Some(r)) = (lay.reference, od.reference) {
            f.fill_rectangle(Point::new(0.0, y), Size::new(lay.seq_x0 - 1.0, h), pal.reference_bg);
            let name = truncate(&od.doc.rows[r].name, ((lay.seq_x0 - 30.0) / 7.0) as usize);
            self.label(f, x, y, h, &name, pal.text, 12.5, true);
        }
        f.fill_rectangle(Point::new(0.0, lay.rows_y0 - 2.0), Size::new(lay.w, 1.0), pal.grid);

        let font_size = (vp.row_h - 3.0).clamp(8.0, 13.0);
        let show_text = vp.row_h >= 8.0;
        let max_chars = ((lay.seq_x0 - 24.0) / (font_size * 0.56)) as usize;
        for (i, y, h) in self.visible_items() {
            match &od.items[i] {
                Item::Group { label, sequences, folded } => {
                    f.fill_rectangle(Point::new(0.0, y), Size::new(lay.seq_x0 - 1.0, h), pal.group_bg);
                    f.fill_rectangle(Point::new(0.0, y), Size::new(lay.w, 1.0), pal.grid);
                    let s = format!("{} {}  ({})", if *folded { "▸" } else { "▾" }, label, group_digits(*sequences));
                    self.label(f, 6.0, y, h, &truncate(&s, ((lay.seq_x0 - 24.0) / 7.0) as usize), pal.text, 12.0, true);
                }
                Item::Seq { row, count } => {
                    if od.selected_rows.contains(row) {
                        f.fill_rectangle(Point::new(0.0, y), Size::new(lay.seq_x0 - 1.0, vp.row_h), pal.row_selected);
                    }
                    if show_text {
                        let indent = if od.display.group_by.is_some() { 18.0 } else { x };
                        let badge = if *count > 1 { format!("×{}", group_digits(*count)) } else { String::new() };
                        let room = max_chars.saturating_sub(badge.len() + 2);
                        let name = truncate(&od.doc.rows[*row].name, room);
                        self.label(f, indent, y, vp.row_h, &name, pal.text, font_size, false);
                        if !badge.is_empty() {
                            f.fill_text(Text {
                                content: badge,
                                position: Point::new(lay.seq_x0 - 8.0, y + vp.row_h / 2.0),
                                color: pal.selection_edge,
                                size: Pixels(font_size - 1.0),
                                align_x: iced::widget::text::Alignment::Right,
                                align_y: Vertical::Center,
                                font: Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT },
                                ..Text::default()
                            });
                        }
                    }
                }
            }
        }
    }

    fn scrollbars(&self, f: &mut Frame) {
        let lay = self.lay;
        let pal = self.pal;
        f.fill_rectangle(Point::new(lay.w - SB, lay.rows_y0), Size::new(SB, lay.rows_h), pal.scrollbar);
        let (ty, th) = self.v.vthumb(lay);
        f.fill(&Path::rounded_rectangle(Point::new(lay.w - SB + 2.0, ty), Size::new(SB - 4.0, th), 4.0.into()), pal.scrollbar_thumb);
        f.fill_rectangle(Point::new(0.0, lay.h - SB), Size::new(lay.w, SB), pal.scrollbar);
        let (tx, tw) = self.v.hthumb(lay);
        f.fill(&Path::rounded_rectangle(Point::new(tx, lay.h - SB + 2.0), Size::new(tw, SB - 4.0), 4.0.into()), pal.scrollbar_thumb);
    }
}

/// Bar height (0..1): linear for identity; for conservation the mismatch
/// share on a log scale (100% and 99.9% full, 99% 2/3, 90% 1/3).
fn graph_height(mode: super::state::GraphMode, v: f32) -> f32 {
    match mode {
        super::state::GraphMode::Identity => v,
        super::state::GraphMode::Conservation => {
            let m = (1.0 - v).max(1e-3);
            (m.log10() / (1e-3f32).log10()).clamp(0.04, 1.0)
        }
    }
}

/// Summary of a pixel bin in the overview strip.
fn bin_value(values: &[Option<f32>], mode: super::state::GraphMode) -> Option<f32> {
    match mode {
        super::state::GraphMode::Identity => values.iter().flatten().copied().reduce(f32::min),
        super::state::GraphMode::Conservation => {
            let (sum, n) = values.iter().flatten().fold((0.0f32, 0u32), |(s, n), v| (s + v, n + 1));
            (n > 0).then(|| sum / n as f32)
        }
    }
}

fn h_off(consensus: bool) -> f32 {
    if consensus { 0.0 } else { 16.0 }
}

/// A 1-2-5 step of at least `min` columns.
pub fn nice_step(min: f32) -> usize {
    let min = min.max(1.0);
    let mut p = 1usize;
    loop {
        for m in [1usize, 2, 5] {
            if (p * m) as f32 >= min {
                return p * m;
            }
        }
        p *= 10;
    }
}

pub fn group_digits(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

fn truncate(s: &str, max: usize) -> String {
    if s.chars().count() <= max {
        s.to_string()
    } else {
        let mut t: String = s.chars().take(max.saturating_sub(1)).collect();
        t.push('…');
        t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers() {
        assert_eq!(nice_step(0.3), 1);
        assert_eq!(nice_step(7.0), 10);
        assert_eq!(nice_step(130.0), 200);
        assert_eq!(group_digits(1234567), "1,234,567");
        assert_eq!(group_digits(12), "12");
        assert_eq!(truncate("abcdef", 4), "abc…");
    }

    #[test]
    fn lanes() {
        let a = |s, e| Annotation { id: 0, name: String::new(), kind: crate::model::AnnotationKind::Region, start: s, end: e, note: String::new() };
        let l = annotation_lanes(&[a(0, 10), a(5, 15), a(10, 20)]);
        assert_eq!(l, vec![0, 1, 0]);
    }
}
