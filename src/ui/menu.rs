//! Menu bar with drop-down menus (File, Edit, Tools).
//!
//! The bar has fixed geometry so the open menu's panel can be pinned right
//! below its title; the panel is drawn as an overlay layer by `view.rs`.

use super::{Message, style};
use iced::widget::{Space, button, column, container, mouse_area, row, rule, text};
use iced::{Alignment, Background, Border, Color, Element, Fill, Font, Length, Shadow, Theme, Vector};

/// Height of the menu bar (also where drop-down panels start).
pub const BAR_H: f32 = 44.0;
const PAD_X: f32 = 12.0;
const BRAND_W: f32 = 132.0;
const TITLE_W: f32 = 64.0;
const SPACING: f32 = 2.0;
const PANEL_W: f32 = 290.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MenuId {
    File,
    Edit,
    Tools,
}

impl MenuId {
    pub const ALL: [MenuId; 3] = [MenuId::File, MenuId::Edit, MenuId::Tools];

    fn label(self) -> &'static str {
        match self {
            MenuId::File => "File",
            MenuId::Edit => "Edit",
            MenuId::Tools => "Tools",
        }
    }

    fn index(self) -> usize {
        MenuId::ALL.iter().position(|&m| m == self).unwrap()
    }
}

/// Left edge of a menu's drop-down panel (under its title).
pub fn title_x(id: MenuId) -> f32 {
    PAD_X + BRAND_W + id.index() as f32 * (TITLE_W + SPACING)
}

pub enum Entry {
    /// `action == None` shows the entry disabled.
    Item { label: &'static str, shortcut: &'static str, action: Option<Message>, danger: bool },
    Separator,
}

impl Entry {
    pub fn item(label: &'static str, shortcut: &'static str, action: Option<Message>) -> Entry {
        Entry::Item { label, shortcut, action, danger: false }
    }

    /// A destructive entry (shown in red).
    pub fn danger(label: &'static str, action: Option<Message>) -> Entry {
        Entry::Item { label, shortcut: "", action, danger: true }
    }

    fn height(&self) -> f32 {
        match self {
            Entry::Item { .. } => 31.0,
            Entry::Separator => 9.0,
        }
    }
}

/// The menu bar; `open` highlights the open menu's title.
pub fn bar<'a>(open: Option<MenuId>, document: Option<&'a str>) -> Element<'a, Message> {
    let mut r = row![
        container(text("PCR Studio").size(17).font(Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT })).width(BRAND_W)
    ]
    .spacing(SPACING)
    .height(BAR_H - 8.0)
    .align_y(Alignment::Center);
    for id in MenuId::ALL {
        let active = open == Some(id);
        let title = button(text(id.label()).size(13).width(Fill).align_x(Alignment::Center))
            .width(TITLE_W)
            .padding([6, 4])
            .style(move |t: &Theme, s| title_style(t, s, active))
            .on_press(Message::ToggleMenu(id));
        // While a menu is open, hovering another title switches to it.
        r = r.push(mouse_area(title).on_enter(Message::HoverMenu(id)));
    }
    r = r.push(Space::new().width(Fill));
    if let Some(name) = document {
        r = r.push(text(name).size(12).style(style::muted));
    }
    container(r).padding([4.0, PAD_X]).width(Fill).height(BAR_H).style(style::toolbar).into()
}

/// The drop-down panel with its entries.
pub fn panel<'a>(entries: Vec<Entry>) -> Element<'a, Message> {
    panel_with_width(entries, PANEL_W)
}

/// Approximate panel height (for keeping context menus inside the window).
pub fn panel_height(entries: &[Entry]) -> f32 {
    entries.iter().map(Entry::height).sum::<f32>() + 10.0
}

pub const CONTEXT_W: f32 = 210.0;

pub fn panel_with_width<'a>(entries: Vec<Entry>, width: f32) -> Element<'a, Message> {
    let mut col = column![].spacing(0);
    for e in entries {
        match e {
            Entry::Separator => col = col.push(container(rule::horizontal(1)).padding([4, 6])),
            Entry::Item { label, shortcut, action, danger } => {
                let enabled = action.is_some();
                let content = row![
                    text(label).size(13).width(Fill),
                    text(shortcut).size(12).style(move |t: &Theme| {
                        let c = t.extended_palette().background.base.text;
                        text::Style { color: Some(Color { a: if enabled { 0.55 } else { 0.3 }, ..c }) }
                    }),
                ]
                .align_y(Alignment::Center);
                col = col.push(
                    button(content)
                        .width(Fill)
                        .padding([6, 12])
                        .style(move |t: &Theme, s| item_style(t, s, danger))
                        .on_press_maybe(action.map(|m| Message::MenuAction(Box::new(m)))),
                );
            }
        }
    }
    container(col).padding(4).width(Length::Fixed(width)).style(panel_style).into()
}

fn title_style(theme: &Theme, status: button::Status, active: bool) -> button::Style {
    let p = theme.extended_palette();
    let bg = if active {
        Some(Background::Color(p.primary.weak.color))
    } else {
        match status {
            button::Status::Hovered => Some(Background::Color(p.background.weak.color)),
            button::Status::Pressed => Some(Background::Color(p.background.strong.color)),
            _ => None,
        }
    };
    button::Style {
        background: bg,
        text_color: if active { p.primary.weak.text } else { p.background.base.text },
        border: Border { color: Color::TRANSPARENT, width: 0.0, radius: 5.0.into() },
        ..Default::default()
    }
}

fn item_style(theme: &Theme, status: button::Status, danger: bool) -> button::Style {
    let p = theme.extended_palette();
    let accent = if danger { p.danger.base } else { p.primary.base };
    let normal = if danger { p.danger.base.color } else { p.background.base.text };
    let (bg, fg) = match status {
        button::Status::Hovered | button::Status::Pressed => (Some(Background::Color(accent.color)), accent.text),
        button::Status::Disabled => (None, Color { a: 0.4, ..normal }),
        _ => (None, normal),
    };
    button::Style { background: bg, text_color: fg, border: Border { radius: 4.0.into(), ..Default::default() }, ..Default::default() }
}

fn panel_style(theme: &Theme) -> container::Style {
    let p = theme.extended_palette();
    container::Style {
        background: Some(Background::Color(p.background.base.color)),
        border: Border { color: p.background.strong.color, width: 1.0, radius: 6.0.into() },
        shadow: Shadow { color: Color { a: 0.22, ..Color::BLACK }, offset: Vector::new(0.0, 4.0), blur_radius: 14.0 },
        text_color: Some(p.background.base.text),
        snap: false,
    }
}
