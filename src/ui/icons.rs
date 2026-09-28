//! Library icons (inline SVG, tinted per kind via the svg style color).

use crate::model::DocKind;
use iced::widget::svg::{self, Handle, Svg};
use iced::{Color, Element, Theme};
use std::sync::LazyLock;

const FOLDER: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="M2.5 6.5A2 2 0 0 1 4.5 4.5h4.4l2.1 2.1h8.5a2 2 0 0 1 2 2v9.4a2 2 0 0 1-2 2h-15a2 2 0 0 1-2-2z" fill="#000"/></svg>"##;

const FOLDER_OPEN: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="M2.5 6.5A2 2 0 0 1 4.5 4.5h4.4l2.1 2.1h7.5a2 2 0 0 1 2 2v.9H7.6a2 2 0 0 0-1.9 1.4L2.5 17.5z" fill="#000" opacity="0.55"/><path d="M5.9 10.9A2 2 0 0 1 7.8 9.5h13.4a1 1 0 0 1 .96 1.28l-2.3 7.9a2 2 0 0 1-1.92 1.42H3.4a1 1 0 0 1-.96-1.28z" fill="#000"/></svg>"##;

/// Sequence list: stacked sequences of different lengths.
const SEQUENCES: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><rect x="3" y="4.5" width="18" height="3.2" rx="1.6" fill="#000"/><rect x="3" y="10.4" width="13" height="3.2" rx="1.6" fill="#000"/><rect x="3" y="16.3" width="16" height="3.2" rx="1.6" fill="#000"/></svg>"##;

/// Alignment: rows of blocks sharing columns.
const ALIGNMENT: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><g fill="#000"><rect x="3" y="4" width="4" height="4" rx="1"/><rect x="8.5" y="4" width="4" height="4" rx="1"/><rect x="14" y="4" width="7" height="4" rx="1"/><rect x="3" y="10" width="4" height="4" rx="1"/><rect x="8.5" y="10" width="4" height="4" rx="1" opacity="0.4"/><rect x="14" y="10" width="7" height="4" rx="1"/><rect x="3" y="16" width="9.5" height="4" rx="1"/><rect x="14" y="16" width="7" height="4" rx="1" opacity="0.4"/></g></svg>"##;

/// Contig: reference on top, mapped reads staggered below.
const CONTIG: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><g fill="#000"><rect x="2" y="3.5" width="20" height="3.6" rx="1.8"/><rect x="3" y="9.6" width="11" height="2.8" rx="1.4" opacity="0.7"/><rect x="9" y="14.1" width="12" height="2.8" rx="1.4" opacity="0.7"/><rect x="5" y="18.6" width="10" height="2.8" rx="1.4" opacity="0.7"/></g></svg>"##;

const CHEVRON_RIGHT: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="M9 5.5l7 6.5-7 6.5" fill="none" stroke="#000" stroke-width="2.6" stroke-linecap="round" stroke-linejoin="round"/></svg>"##;

const CHEVRON_DOWN: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="M5.5 9l6.5 7 6.5-7" fill="none" stroke="#000" stroke-width="2.6" stroke-linecap="round" stroke-linejoin="round"/></svg>"##;

fn handle(src: &'static str) -> Handle {
    Handle::from_memory(src.as_bytes())
}

static FOLDER_H: LazyLock<Handle> = LazyLock::new(|| handle(FOLDER));
static FOLDER_OPEN_H: LazyLock<Handle> = LazyLock::new(|| handle(FOLDER_OPEN));
static SEQUENCES_H: LazyLock<Handle> = LazyLock::new(|| handle(SEQUENCES));
static ALIGNMENT_H: LazyLock<Handle> = LazyLock::new(|| handle(ALIGNMENT));
static CONTIG_H: LazyLock<Handle> = LazyLock::new(|| handle(CONTIG));
static CHEVRON_RIGHT_H: LazyLock<Handle> = LazyLock::new(|| handle(CHEVRON_RIGHT));
static CHEVRON_DOWN_H: LazyLock<Handle> = LazyLock::new(|| handle(CHEVRON_DOWN));

/// Expand/collapse chevron in the theme's muted text color.
pub fn chevron<'a, M: 'a>(expanded: bool) -> Element<'a, M> {
    Svg::new(if expanded { CHEVRON_DOWN_H.clone() } else { CHEVRON_RIGHT_H.clone() })
        .width(12)
        .height(12)
        .style(|t: &Theme, _s| svg::Style { color: Some(Color { a: 0.7, ..t.extended_palette().background.base.text }) })
        .into()
}

fn icon<'a, M: 'a>(h: &Handle, size: f32, color: Color) -> Element<'a, M> {
    Svg::new(h.clone()).width(size).height(size).style(move |_t: &Theme, _s| svg::Style { color: Some(color) }).into()
}

pub fn folder_color(dark: bool) -> Color {
    if dark { Color::from_rgb8(222, 170, 70) } else { Color::from_rgb8(221, 157, 36) }
}

pub fn folder<'a, M: 'a>(open: bool, dark: bool) -> Element<'a, M> {
    icon(if open { &FOLDER_OPEN_H } else { &FOLDER_H }, 18.0, folder_color(dark))
}

pub fn document<'a, M: 'a>(kind: DocKind, dark: bool) -> Element<'a, M> {
    let (h, light, darkc) = match kind {
        DocKind::Sequences => (&*SEQUENCES_H, Color::from_rgb8(58, 118, 214), Color::from_rgb8(110, 160, 240)),
        DocKind::Alignment => (&*ALIGNMENT_H, Color::from_rgb8(22, 150, 136), Color::from_rgb8(70, 190, 170)),
        DocKind::Contig => (&*CONTIG_H, Color::from_rgb8(140, 82, 204), Color::from_rgb8(180, 130, 240)),
    };
    icon(h, 17.0, if dark { darkc } else { light })
}
