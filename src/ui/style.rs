//! Colors and widget styles.

use iced::widget::{button, container};
use iced::{Background, Border, Color, Shadow, Theme, Vector};

pub const fn rgb(r: u8, g: u8, b: u8) -> Color {
    Color { r: r as f32 / 255.0, g: g as f32 / 255.0, b: b as f32 / 255.0, a: 1.0 }
}

pub const fn rgba(r: u8, g: u8, b: u8, a: f32) -> Color {
    Color { r: r as f32 / 255.0, g: g as f32 / 255.0, b: b as f32 / 255.0, a }
}

/// Colors used by the alignment canvas.
#[derive(Debug, Clone, Copy)]
pub struct Palette {
    pub bg: Color,
    pub names_bg: Color,
    pub header_bg: Color,
    pub text: Color,
    pub text_muted: Color,
    pub grid: Color,
    pub base_a: Color,
    pub base_c: Color,
    pub base_g: Color,
    pub base_t: Color,
    pub base_other: Color,
    pub base_text: Color,
    pub same: Color,
    pub same_text: Color,
    pub gap: Color,
    pub gap_strong: Color,
    pub reference_bg: Color,
    pub consensus_bg: Color,
    pub id_high: Color,
    pub id_mid: Color,
    pub id_low: Color,
    /// Conservation tiers: 100%, >=99%, >=95%, >=90%, >=75%.
    pub cons: [Color; 5],
    pub selection: Color,
    pub selection_edge: Color,
    pub hit: Color,
    pub group_bg: Color,
    pub row_selected: Color,
    pub scrollbar: Color,
    pub scrollbar_thumb: Color,
    pub ann_fwd: Color,
    pub ann_rev: Color,
    pub ann_probe: Color,
    pub ann_region: Color,
    /// Base quality Q0..=Q40: red, neutral gray at Q20, blue ([`quality_scale`]).
    pub quality: [Color; QUALITY_MAX as usize + 1],
}

/// Highest Phred quality with its own color; higher ones look the same.
pub const QUALITY_MAX: u8 = 40;

/// Diverging scale for base qualities: `low` at Q0, `mid` at Q20, `high` from Q40.
fn quality_scale(low: Color, mid: Color, high: Color) -> [Color; QUALITY_MAX as usize + 1] {
    let half = QUALITY_MAX as f32 / 2.0;
    let mix = |a: Color, b: Color, t: f32| Color { r: a.r + (b.r - a.r) * t, g: a.g + (b.g - a.g) * t, b: a.b + (b.b - a.b) * t, a: 1.0 };
    std::array::from_fn(|q| {
        let q = q as f32;
        if q <= half { mix(low, mid, q / half) } else { mix(mid, high, (q - half) / half) }
    })
}

impl Palette {
    pub fn light() -> Palette {
        Palette {
            bg: rgb(255, 255, 255),
            names_bg: rgb(247, 248, 250),
            header_bg: rgb(242, 244, 248),
            text: rgb(28, 32, 40),
            text_muted: rgb(110, 118, 130),
            grid: rgb(222, 226, 232),
            base_a: rgb(236, 104, 98),
            base_c: rgb(92, 142, 230),
            base_g: rgb(242, 183, 64),
            base_t: rgb(92, 178, 104),
            base_other: rgb(178, 182, 190),
            base_text: rgb(20, 22, 28),
            same: rgb(232, 235, 240),
            same_text: rgb(150, 156, 166),
            gap: rgb(170, 175, 184),
            gap_strong: rgb(40, 40, 48),
            reference_bg: rgb(255, 247, 204),
            consensus_bg: rgb(234, 240, 250),
            id_high: rgb(46, 166, 76),
            id_mid: rgb(186, 190, 64),
            id_low: rgb(214, 70, 70),
            cons: [rgb(30, 140, 70), rgb(96, 186, 96), rgb(186, 204, 70), rgb(236, 190, 60), rgb(238, 132, 52)],
            selection: rgba(40, 110, 240, 0.16),
            selection_edge: rgb(40, 110, 240),
            hit: rgb(230, 40, 160),
            group_bg: rgb(228, 233, 242),
            row_selected: rgba(40, 110, 240, 0.14),
            scrollbar: rgb(238, 240, 244),
            scrollbar_thumb: rgb(186, 192, 202),
            ann_fwd: rgb(88, 160, 240),
            ann_rev: rgb(242, 130, 64),
            ann_probe: rgb(160, 104, 220),
            ann_region: rgb(120, 190, 150),
            quality: quality_scale(rgb(190, 40, 45), rgb(226, 225, 221), rgb(28, 92, 171)),
        }
    }

    pub fn dark() -> Palette {
        Palette {
            bg: rgb(28, 30, 36),
            names_bg: rgb(34, 37, 44),
            header_bg: rgb(38, 41, 49),
            text: rgb(226, 230, 236),
            text_muted: rgb(140, 148, 160),
            grid: rgb(58, 62, 72),
            base_a: rgb(200, 84, 80),
            base_c: rgb(74, 120, 206),
            base_g: rgb(210, 158, 52),
            base_t: rgb(74, 156, 88),
            base_other: rgb(110, 116, 126),
            base_text: rgb(245, 246, 248),
            same: rgb(52, 56, 66),
            same_text: rgb(120, 128, 140),
            gap: rgb(100, 106, 118),
            gap_strong: rgb(236, 236, 240),
            reference_bg: rgb(74, 66, 30),
            consensus_bg: rgb(40, 50, 68),
            id_high: rgb(50, 160, 80),
            id_mid: rgb(160, 164, 60),
            id_low: rgb(200, 70, 70),
            cons: [rgb(40, 150, 80), rgb(80, 170, 86), rgb(160, 176, 60), rgb(200, 160, 50), rgb(206, 112, 44)],
            selection: rgba(90, 150, 255, 0.22),
            selection_edge: rgb(110, 165, 255),
            hit: rgb(255, 80, 190),
            group_bg: rgb(46, 52, 64),
            row_selected: rgba(90, 150, 255, 0.20),
            scrollbar: rgb(36, 39, 46),
            scrollbar_thumb: rgb(86, 92, 104),
            ann_fwd: rgb(70, 136, 220),
            ann_rev: rgb(220, 112, 50),
            ann_probe: rgb(140, 90, 200),
            ann_region: rgb(90, 160, 120),
            quality: quality_scale(rgb(230, 103, 103), rgb(72, 72, 70), rgb(109, 167, 236)),
        }
    }

    #[inline]
    pub fn base(&self, b: u8) -> Color {
        match b {
            b'A' => self.base_a,
            b'C' => self.base_c,
            b'G' => self.base_g,
            b'T' => self.base_t,
            _ => self.base_other,
        }
    }

    #[inline]
    pub fn quality(&self, q: u8) -> Color {
        self.quality[q.min(QUALITY_MAX) as usize]
    }

    /// Text on a cell of quality `q`: dark or light, whichever contrasts more.
    pub fn quality_text(&self, q: u8) -> Color {
        if luminance(self.quality(q)) > 0.2 { DARK_TEXT } else { LIGHT_TEXT }
    }

    pub fn identity(&self, id: f32) -> Color {
        if id >= 0.9999 {
            self.id_high
        } else if id >= 0.3 {
            self.id_mid
        } else {
            self.id_low
        }
    }

    pub fn conservation(&self, c: f32) -> Color {
        match c {
            c if c >= 0.9999 => self.cons[0],
            c if c >= 0.99 => self.cons[1],
            c if c >= 0.95 => self.cons[2],
            c if c >= 0.90 => self.cons[3],
            c if c >= 0.75 => self.cons[4],
            _ => self.id_low,
        }
    }

    pub fn graph(&self, mode: super::state::GraphMode, v: f32) -> Color {
        match mode {
            super::state::GraphMode::Conservation => self.conservation(v),
            super::state::GraphMode::Identity => self.identity(v),
        }
    }

    /// Color of an annotation: its own, else the color of its type.
    pub fn annotation(&self, a: &crate::model::Annotation) -> Color {
        a.color.map_or_else(|| self.annotation_kind(a.kind), from_rgb)
    }

    pub fn annotation_kind(&self, kind: crate::model::AnnotationKind) -> Color {
        use crate::model::AnnotationKind::*;
        match kind {
            ForwardPrimer => self.ann_fwd,
            ReversePrimer => self.ann_rev,
            Probe => self.ann_probe,
            Region => self.ann_region,
        }
    }

    pub fn for_theme(theme: &Theme) -> Palette {
        if theme.extended_palette().is_dark { Palette::dark() } else { Palette::light() }
    }
}

const DARK_TEXT: Color = rgb(20, 22, 28);
const LIGHT_TEXT: Color = rgb(245, 246, 248);

/// Relative luminance (WCAG) of an sRGB color.
fn luminance(c: Color) -> f32 {
    let lin = |v: f32| if v <= 0.04045 { v / 12.92 } else { ((v + 0.055) / 1.055).powf(2.4) };
    0.2126 * lin(c.r) + 0.7152 * lin(c.g) + 0.0722 * lin(c.b)
}

/// Label on an annotation of color `c`: white, dark on light colors (the
/// colors of the annotation types all keep white labels).
pub fn annotation_text(c: Color) -> Color {
    if luminance(c) > 0.5 { DARK_TEXT } else { Color::WHITE }
}

pub const fn from_rgb(c: crate::model::Rgb) -> Color {
    rgb(c.0, c.1, c.2)
}

/// Preset colors for annotations, readable in light and dark mode.
pub const ANNOTATION_COLORS: [crate::model::Rgb; 12] = {
    use crate::model::Rgb;
    [
        Rgb(0xD6, 0x3C, 0x3C), // red
        Rgb(0xE8, 0x83, 0x2E), // orange
        Rgb(0xF2, 0xC1, 0x2E), // yellow
        Rgb(0x8B, 0xC3, 0x4A), // lime
        Rgb(0x2E, 0x9E, 0x5B), // green
        Rgb(0x1F, 0xA3, 0x9B), // teal
        Rgb(0x3B, 0xA7, 0xE0), // sky
        Rgb(0x34, 0x5F, 0xD0), // blue
        Rgb(0x7B, 0x5C, 0xD6), // violet
        Rgb(0xD9, 0x4F, 0xA4), // pink
        Rgb(0x8D, 0x62, 0x44), // brown
        Rgb(0x7E, 0x86, 0x92), // gray
    ]
};

pub fn app_theme(dark: bool) -> Theme {
    if dark {
        Theme::custom(
            "PCR Studio Dark",
            iced::theme::Palette {
                background: rgb(30, 32, 38),
                text: rgb(226, 230, 236),
                primary: rgb(76, 136, 240),
                success: rgb(60, 170, 90),
                warning: rgb(230, 170, 50),
                danger: rgb(220, 80, 80),
            },
        )
    } else {
        Theme::custom(
            "PCR Studio",
            iced::theme::Palette {
                background: rgb(250, 251, 253),
                text: rgb(28, 32, 40),
                primary: rgb(38, 104, 220),
                success: rgb(40, 150, 70),
                warning: rgb(220, 150, 30),
                danger: rgb(200, 60, 60),
            },
        )
    }
}

pub fn panel(theme: &Theme) -> container::Style {
    let p = theme.extended_palette();
    container::Style {
        background: Some(Background::Color(p.background.weak.color)),
        border: Border { color: p.background.strong.color, width: 0.0, radius: 0.0.into() },
        ..Default::default()
    }
}

pub fn toolbar(theme: &Theme) -> container::Style {
    let p = theme.extended_palette();
    container::Style {
        background: Some(Background::Color(p.background.base.color)),
        border: Border { color: p.background.strong.color, width: 1.0, radius: 0.0.into() },
        ..Default::default()
    }
}

pub fn card(theme: &Theme) -> container::Style {
    let p = theme.extended_palette();
    container::Style {
        background: Some(Background::Color(p.background.base.color)),
        border: Border { color: p.background.strong.color, width: 1.0, radius: 8.0.into() },
        shadow: Shadow { color: rgba(0, 0, 0, 0.25), offset: Vector::new(0.0, 4.0), blur_radius: 18.0 },
        text_color: Some(p.background.base.text),
        snap: false,
    }
}

pub fn backdrop(_theme: &Theme) -> container::Style {
    container::Style { background: Some(Background::Color(rgba(0, 0, 0, 0.35))), ..Default::default() }
}

pub fn muted(theme: &Theme) -> iced::widget::text::Style {
    let text = theme.extended_palette().background.base.text;
    iced::widget::text::Style { color: Some(Color { a: 0.62, ..text }) }
}

/// Text in the danger color (missing files, failures).
pub fn danger_text(theme: &Theme) -> iced::widget::text::Style {
    iced::widget::text::Style { color: Some(theme.extended_palette().danger.base.color) }
}

/// Text in the success color (checked and fine).
pub fn success_text(theme: &Theme) -> iced::widget::text::Style {
    iced::widget::text::Style { color: Some(theme.extended_palette().success.base.color) }
}

/// Flat toolbar button.
pub fn tool_button(theme: &Theme, status: button::Status) -> button::Style {
    let p = theme.extended_palette();
    let bg = match status {
        button::Status::Hovered => Some(Background::Color(p.background.weak.color)),
        button::Status::Pressed => Some(Background::Color(p.background.strong.color)),
        _ => None,
    };
    button::Style {
        background: bg,
        text_color: if status == button::Status::Disabled { p.background.strong.color } else { p.background.base.text },
        border: Border { color: Color::TRANSPARENT, width: 0.0, radius: 5.0.into() },
        ..Default::default()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::AnnotationKind;

    #[test]
    fn annotation_labels() {
        for pal in [Palette::light(), Palette::dark()] {
            for kind in AnnotationKind::ALL {
                assert_eq!(annotation_text(pal.annotation_kind(kind)), Color::WHITE, "{kind}");
            }
        }
        assert_eq!(annotation_text(from_rgb(ANNOTATION_COLORS[2])), DARK_TEXT);
        assert_eq!(annotation_text(from_rgb(ANNOTATION_COLORS[7])), Color::WHITE);
    }
}
