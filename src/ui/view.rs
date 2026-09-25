//! Layout of the main window.

use super::canvas::{AlignmentView, group_digits};
use super::dialogs::{self, ModalMsg};
use super::state::{BASE_COL_W, GraphMode, Highlight, OpenDoc};
use super::{App, FolderChoice, GroupChoice, Message, style};
use crate::display::{Item, SortKey};
use crate::model::{AnnotationKind, ConsensusThreshold, DocKind};
use crate::primer::{self, OligoSource, Orientation, PrimerReport};
use crate::search::Scope;
use iced::widget::canvas::{self as cv, Frame, Geometry};
use iced::widget::{
    Space, button, canvas, center, checkbox, column, container, mouse_area, opaque, pick_list, progress_bar, row, rule, scrollable, slider,
    stack, text, text_input, tooltip,
};
use iced::{Alignment, Color, Element, Fill, Font, Length, Point, Rectangle, Renderer, Size, Theme, mouse};
use std::collections::HashSet;

const LIBRARY_W: f32 = 260.0;
const OPTIONS_W: f32 = 300.0;

fn tool<'a>(label: &'a str, msg: Option<Message>, tip: &'a str) -> Element<'a, Message> {
    let b = button(text(label).size(13)).padding([6, 10]).style(style::tool_button).on_press_maybe(msg);
    tooltip(b, container(text(tip).size(12)).padding(6).style(container::rounded_box), tooltip::Position::Bottom).into()
}

fn section<'a>(title: &'a str) -> Element<'a, Message> {
    column![text(title).size(13).font(Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT }), rule::horizontal(1)].spacing(4).into()
}

impl App {
    pub(super) fn view(&self) -> Element<'_, Message> {
        let body = row![self.library_panel(), rule::vertical(1), self.center_panel()];
        let body: Element<Message> = if self.open.is_some() { row![body, rule::vertical(1), self.options_panel()].into() } else { body.into() };
        let main = column![self.toolbar(), container(body).height(Fill), self.status_bar()];
        let base: Element<Message> = container(main).width(Fill).height(Fill).into();
        if let Some(m) = &self.modal {
            let overlay = opaque(mouse_area(center(opaque(dialogs::view_modal(m))).style(style::backdrop)).on_press(Message::Noop));
            stack![base, overlay].into()
        } else {
            base
        }
    }

    fn toolbar(&self) -> Element<'_, Message> {
        let idle = !self.busy();
        let has_doc = self.open.is_some();
        let bar = row![
            text("PCR Studio").size(17).font(Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT }),
            Space::new().width(18),
            tool("Import FASTA", idle.then_some(Message::ImportFasta), "Import one or more FASTA files into the library (Ctrl+O)"),
            tool("Export", (idle && has_doc).then_some(Message::OpenExport), "Export sequences, selection or consensus as FASTA (Ctrl+E)"),
            rule::vertical(1),
            tool("Map to Reference", (idle && has_doc).then_some(Message::OpenMapDialog), "Map the sequences of this document to a reference (Ctrl+M)"),
            rule::vertical(1),
            tool("Metadata table", (idle && has_doc).then_some(Message::ImportMetadata), "Import a CSV/TSV table: first column = sequence name"),
            tool("Names → fields", (idle && has_doc).then_some(Message::OpenNameFields), "Split sequence names into metadata fields"),
            Space::new().width(Fill),
            tool("Demo data", idle.then_some(Message::DemoData), "Add a synthetic 35 kb × 2,000 genome data set to the library"),
            tool("Settings", Some(Message::OpenSettings), "Tm conditions, 3' window, theme"),
        ]
        .spacing(4)
        .height(36)
        .align_y(Alignment::Center);
        container(bar).padding([4, 12]).width(Fill).style(style::toolbar).into()
    }

    fn library_panel(&self) -> Element<'_, Message> {
        let known: HashSet<i64> = self.folders.iter().map(|f| f.id).collect();
        let mut entries: Vec<Element<Message>> = Vec::new();
        self.push_tree(&mut entries, None, 0, &known);
        let mut list = column(entries).spacing(2);
        if self.library.is_empty() && self.folders.is_empty() {
            list = list.push(text("No documents yet. Import FASTA files or generate demo data.").size(12).style(style::muted));
        }
        let idle = !self.busy();
        let (open, rename, delete) = match (self.library_selected, self.folder_selected) {
            (Some(d), _) => (idle.then_some(Message::OpenDocument(d)), Some(Message::AskRename(d)), idle.then_some(Message::AskDelete(d))),
            (None, Some(f)) => (None, Some(Message::AskRenameFolder(f)), Some(Message::AskDeleteFolder(f))),
            _ => (None, None, None),
        };
        let actions = row![
            button(text("Open").size(12)).padding([4, 10]).on_press_maybe(open),
            button(text("Rename").size(12)).padding([4, 10]).style(button::secondary).on_press_maybe(rename),
            button(text("Delete").size(12)).padding([4, 10]).style(button::danger).on_press_maybe(delete),
        ]
        .spacing(6);
        let mut organize = row![
            tooltip(
                button(text("New folder").size(12)).padding([4, 10]).style(button::secondary).on_press(Message::NewFolder),
                container(text("Creates a folder inside the selected folder (or next to the selected document)").size(12))
                    .padding(6)
                    .style(container::rounded_box),
                tooltip::Position::Top,
            ),
        ]
        .spacing(6)
        .align_y(Alignment::Center);
        if self.library_selected.is_some() || self.folder_selected.is_some() {
            organize = organize.push(
                pick_list(self.move_targets(), None::<FolderChoice>, Message::MoveTo).placeholder("Move to…").text_size(12).padding(4).width(Fill),
            );
        }
        let info: Option<Element<Message>> = if let Some(d) = self.library_selected.and_then(|id| self.library.iter().find(|d| d.id == id)) {
            Some(
                column![
                    text(d.description.as_str()).size(11).style(style::muted).wrapping(text::Wrapping::WordOrGlyph),
                    text(format!("Created {}", d.created)).size(11).style(style::muted)
                ]
                .spacing(2)
                .into(),
            )
        } else {
            self.folder_selected.map(|id| {
                let n = self.folder_doc_count(id);
                text(format!("Folder · {n} document(s) · click again to expand/collapse")).size(11).style(style::muted).into()
            })
        };
        let mut col = column![
            row![text("Library").size(15), Space::new().width(Fill), text(format!("{}", self.library.len())).size(12).style(style::muted)],
            scrollable(list).height(Fill),
            actions,
            organize,
        ]
        .spacing(8);
        if let Some(info) = info {
            col = col.push(info);
        }
        container(col).padding(10).width(LIBRARY_W).height(Fill).style(style::panel).into()
    }

    /// Adds the folders and documents below `parent` (depth-first) to `out`.
    fn push_tree<'a>(&'a self, out: &mut Vec<Element<'a, Message>>, parent: Option<i64>, depth: usize, known: &HashSet<i64>) {
        let indent = depth as f32 * 14.0;
        // Unknown parents (should not happen) fall back to the top level.
        let parent_of = |p: Option<i64>| p.filter(|id| known.contains(id));
        for f in self.folders.iter().filter(|f| parent_of(f.parent) == parent) {
            let selected = self.folder_selected == Some(f.id);
            let arrow = button(text(if f.collapsed { "▸" } else { "▾" }).size(18).line_height(1.0))
                .padding([2, 4])
                .style(button::text)
                .on_press(Message::ToggleFolder(f.id));
            let label = row![
                text(f.name.as_str()).size(13).font(Font { weight: iced::font::Weight::Semibold, ..Font::DEFAULT }).width(Fill),
                text(format!("{}", self.folder_doc_count(f.id))).size(11).style(style::muted),
            ]
            .align_y(Alignment::Center);
            out.push(
                row![
                    Space::new().width(indent),
                    arrow,
                    button(label).width(Fill).padding([5, 6]).style(style::list_item(selected)).on_press(Message::SelectFolder(f.id)),
                ]
                .align_y(Alignment::Center)
                .into(),
            );
            if !f.collapsed {
                self.push_tree(out, Some(f.id), depth + 1, known);
            }
        }
        for d in self.library.iter().filter(|d| parent_of(d.folder) == parent) {
            let selected = self.library_selected == Some(d.id);
            let open = self.open.as_ref().is_some_and(|o| o.info.id == d.id);
            let kind = match d.kind {
                DocKind::Contig => "Contig",
                DocKind::Alignment => "Alignment",
                DocKind::Sequences => "Sequences",
            };
            let subtitle = format!("{kind} · {} seqs · {} bp", group_digits(d.n_rows), group_digits(d.width));
            let name = if open { format!("● {}", d.name) } else { d.name.clone() };
            let entry = column![text(name).size(13), text(subtitle).size(11).style(style::muted)].spacing(1);
            // Documents line up with folder names (past the arrow button).
            let pad = if depth > 0 || !self.folders.is_empty() { indent + 20.0 } else { indent };
            out.push(
                row![
                    Space::new().width(pad),
                    button(entry).width(Fill).padding([5, 8]).style(style::list_item(selected)).on_press(Message::SelectLibraryDoc(d.id)),
                ]
                .into(),
            );
        }
    }

    /// Documents in a folder including its subfolders.
    fn folder_doc_count(&self, id: i64) -> usize {
        let sub: HashSet<i64> = self.folder_subtree(id).into_iter().collect();
        self.library.iter().filter(|d| d.folder.is_some_and(|f| sub.contains(&f))).count()
    }

    /// Destinations for "Move to": top level plus all folders in tree order,
    /// without the selected folder's own subtree.
    fn move_targets(&self) -> Vec<FolderChoice> {
        let exclude: HashSet<i64> = self.folder_selected.map(|f| self.folder_subtree(f)).unwrap_or_default().into_iter().collect();
        let mut out = vec![FolderChoice { id: None, label: "Library (top level)".into() }];
        fn walk(app: &App, parent: Option<i64>, depth: usize, exclude: &HashSet<i64>, out: &mut Vec<FolderChoice>) {
            for f in app.folders.iter().filter(|f| f.parent == parent && !exclude.contains(&f.id)) {
                out.push(FolderChoice { id: Some(f.id), label: format!("{}{}", "· ".repeat(depth + 1), f.name) });
                walk(app, Some(f.id), depth + 1, exclude, out);
            }
        }
        walk(self, None, 0, &exclude, &mut out);
        out
    }

    fn center_panel(&self) -> Element<'_, Message> {
        let Some(od) = &self.open else {
            return self.welcome();
        };
        let viewer = canvas(AlignmentView { od, prefs: &self.settings.view, pal: self.palette() }).width(Fill).height(Fill);
        let mut layers: Vec<Element<Message>> = vec![viewer.into()];
        if od.popup_open
            && let Some(r) = &od.report {
                layers.push(
                    container(self.primer_popup(od, r))
                        .width(Fill)
                        .height(Fill)
                        .align_x(Alignment::End)
                        .align_y(Alignment::Start)
                        .padding(iced::Padding { top: 64.0, right: 22.0, bottom: 20.0, left: 0.0 })
                        .into(),
                );
            }
        column![self.viewer_toolbar(od), iced::widget::Stack::with_children(layers).width(Fill).height(Fill)].into()
    }

    fn welcome(&self) -> Element<'_, Message> {
        let content = column![
            text("PCR Studio").size(28),
            text("Sequence alignment viewer for PCR design and inclusivity evaluation").size(15).style(style::muted),
            Space::new().height(12),
            text("1. Import FASTA files (aligned or unaligned) — they are stored in the library.").size(14),
            text("2. Open a document from the library (click twice or press Open).").size(14),
            text("3. Map sequences to a reference, or view an existing alignment.").size(14),
            text("4. Drag across columns to evaluate a primer or probe against all sequences.").size(14),
            Space::new().height(12),
            row![
                button("Import FASTA…").on_press_maybe((!self.busy()).then_some(Message::ImportFasta)).style(button::primary),
                button("Generate demo data").on_press_maybe((!self.busy()).then_some(Message::DemoData)).style(button::secondary),
            ]
            .spacing(10),
        ]
        .spacing(6)
        .max_width(620);
        center(content).into()
    }

    fn viewer_toolbar<'a>(&'a self, od: &'a OpenDoc) -> Element<'a, Message> {
        let zoom_pct = od.vp.col_w / BASE_COL_W * 100.0;
        let zoom_label = if zoom_pct >= 10.0 { format!("{zoom_pct:.0}%") } else { format!("{zoom_pct:.2}%") };
        let scopes: Vec<Scope> = if od.reference.is_some() {
            vec![Scope::Reference, Scope::Consensus, Scope::AllSequences]
        } else {
            vec![Scope::Consensus, Scope::AllSequences]
        };
        let hits = &od.search;
        let hit_label = match (hits.current, hits.hits.len()) {
            (_, 0) if hits.searched => "0 hits".to_string(),
            (Some(c), n) if n > 0 => format!("{} / {}", group_digits(c + 1), group_digits(n)),
            _ => String::new(),
        };
        let jump_hint = if od.reference.is_some() { "Ref position or range" } else { "Column or range" };
        let bar = row![
            tool("−", Some(Message::ZoomOut), "Zoom out (Ctrl+wheel, Ctrl+-)"),
            text(zoom_label).size(12).width(52).align_x(Alignment::Center),
            tool("+", Some(Message::ZoomIn), "Zoom in (Ctrl+wheel, Ctrl+=)"),
            tool("Fit", Some(Message::ZoomFit), "Show the whole alignment width (Ctrl+0)"),
            tool("1:1", Some(Message::ZoomLetters), "Zoom to readable bases"),
            rule::vertical(1),
            tooltip(
                row![text("Rows").size(12), slider(2.0..=24.0, od.vp.row_h, Message::RowHeight).step(1.0).width(70)]
                    .spacing(6)
                    .align_y(Alignment::Center),
                container(text("Row height: small values show thousands of sequences at once").size(12)).padding(6).style(container::rounded_box),
                tooltip::Position::Bottom,
            ),
            rule::vertical(1),
            text_input(jump_hint, &od.jump_text).on_input(Message::JumpInput).on_submit(Message::JumpSubmit).padding(5).size(13).width(130),
            tool("Go", Some(Message::JumpSubmit), "Jump to position; a range (e.g. 100-120) is selected"),
            rule::vertical(1),
            text_input("Search (IUPAC)…", &od.search.query)
                .on_input(Message::SearchInput)
                .on_submit(Message::SearchSubmit)
                .padding(5)
                .size(13)
                .font(Font::MONOSPACE)
                .width(Length::FillPortion(1)),
            pick_list(scopes, Some(od.search.scope), Message::SearchScope).text_size(12).padding(5),
            pick_list(
                (0u8..=3).map(Mismatches).collect::<Vec<_>>(),
                Some(Mismatches(od.search.mismatches)),
                |m: Mismatches| Message::SearchMismatches(m.0)
            )
            .text_size(12)
            .padding(5)
            .width(70),
            tool("‹", Some(Message::SearchStep(-1)), "Previous hit"),
            tool("›", Some(Message::SearchStep(1)), "Next hit (Enter)"),
            text(hit_label).size(12).width(Length::Shrink),
        ]
        .spacing(6)
        .align_y(Alignment::Center)
        .height(38);
        container(bar).padding([2, 10]).width(Fill).style(style::toolbar).into()
    }

    fn options_panel(&self) -> Element<'_, Message> {
        let od = self.open.as_ref().unwrap();
        let prefs = &self.settings.view;
        let doc = &od.doc;
        let kind = od.info.kind.label();
        let ref_name = od.reference.map(|r| doc.rows[r].name.clone()).unwrap_or_else(|| "none".into());
        let one_selected = (od.selected_rows.len() == 1).then(|| *od.selected_rows.iter().next().unwrap());

        let doc_info = column![
            section("Document"),
            text(od.info.name.as_str()).size(14),
            text(format!("{kind} · {} sequences · {} columns", group_digits(doc.rows.len()), group_digits(doc.width))).size(12).style(style::muted),
            text(format!("Reference: {ref_name}")).size(12),
            row![
                button(text("Set selected row as reference").size(12))
                    .padding([4, 8])
                    .style(button::secondary)
                    .on_press_maybe(one_selected.map(|_| Message::SetReferenceFromSelection)),
                button(text("Clear").size(12)).padding([4, 8]).style(button::text).on_press_maybe(od.reference.map(|_| Message::ClearReference)),
            ]
            .spacing(6),
            if od.selected_rows.is_empty() {
                Element::from(text("Click names to select rows (Ctrl/Shift for several).").size(11).style(style::muted))
            } else {
                row![
                    text(format!("{} row(s) selected", od.selected_rows.len())).size(12),
                    button(text("Clear").size(11)).padding([2, 6]).style(button::text).on_press(Message::ClearRowSelection)
                ]
                .spacing(6)
                .align_y(Alignment::Center)
                .into()
            },
        ]
        .spacing(6);

        let display = column![
            section("Display"),
            checkbox(prefs.show_consensus).label("Consensus").on_toggle(Message::ToggleConsensus).text_size(13),
            row![
                text("Threshold").size(12).width(80),
                pick_list(&ConsensusThreshold::ALL[..], Some(prefs.threshold), Message::SetThreshold).text_size(12).padding(4).width(Fill)
            ]
            .align_y(Alignment::Center),
            checkbox(prefs.show_identity).label("Graph").on_toggle(Message::ToggleIdentity).text_size(13),
            pick_list([GraphMode::Conservation, GraphMode::Identity], Some(prefs.graph), Message::SetGraph).text_size(12).padding(4).width(Fill),
            checkbox(prefs.show_overview).label("Overview strip").on_toggle(Message::ToggleOverview).text_size(13),
            checkbox(prefs.show_annotations).label("Annotations").on_toggle(Message::ToggleAnnotations).text_size(13),
            checkbox(prefs.color_bases).label("Color bases").on_toggle(Message::ToggleColors).text_size(13),
        ]
        .spacing(6);

        let mut modes = vec![Highlight::None, Highlight::Consensus];
        if od.reference.is_some() {
            modes.insert(1, Highlight::Reference);
        }
        let effective = if prefs.highlight == Highlight::Reference && od.reference.is_none() { Highlight::None } else { prefs.highlight };
        let highlighting = column![
            section("Highlighting"),
            pick_list(modes, Some(effective), Message::SetHighlight).text_size(12).padding(4).width(Fill),
            checkbox(prefs.use_dots).label("Use dots for identical bases").on_toggle(Message::ToggleDots).text_size(13),
            checkbox(prefs.highlight_gaps).label("Highlight gaps").on_toggle(Message::ToggleGaps).text_size(13),
        ]
        .spacing(6);

        let mut sort_keys = vec![SortKey::Original, SortKey::Name, SortKey::StartPosition, SortKey::Length];
        if od.reference.is_some() {
            sort_keys.push(SortKey::IdentityToReference);
        }
        sort_keys.extend(od.meta_keys.iter().map(|k| SortKey::Meta(k.clone())));
        let mut groups = vec![GroupChoice::None];
        groups.extend(od.meta_keys.iter().map(|k| GroupChoice::Key(k.clone())));
        let group = match &od.display.group_by {
            None => GroupChoice::None,
            Some(k) => GroupChoice::Key(k.clone()),
        };
        let shown = od.items.iter().filter(|i| matches!(i, Item::Seq { .. })).count();
        let rows = column![
            section("Sequences"),
            text_input("Filter by name…", &od.display.name_filter).on_input(Message::FilterInput).padding(5).size(13),
            row![text("Sort").size(12).width(60), pick_list(sort_keys, Some(od.display.sort.clone()), Message::SortBy).text_size(12).padding(4).width(Fill)]
                .align_y(Alignment::Center),
            checkbox(od.display.descending).label("Descending").on_toggle(Message::SortDescending).text_size(13),
            row![text("Group").size(12).width(60), pick_list(groups, Some(group), Message::GroupBy).text_size(12).padding(4).width(Fill)]
                .align_y(Alignment::Center),
            checkbox(od.display.collapse_identical).label("Collapse identical sequences").on_toggle(Message::Collapse).text_size(13),
            text(format!("{} rows shown", group_digits(shown))).size(11).style(style::muted),
            if od.meta_keys.is_empty() {
                Element::from(text("No metadata yet: import a table or split names into fields (toolbar).").size(11).style(style::muted))
            } else {
                text(format!("Metadata: {}", od.meta_keys.join(", "))).size(11).style(style::muted).into()
            },
        ]
        .spacing(6);

        let mut anns = column![section("Annotations")].spacing(4);
        if od.annotations.is_empty() {
            anns = anns.push(text("Select columns and use “Add annotation” in the oligo report.").size(11).style(style::muted));
        }
        for (i, a) in od.annotations.iter().enumerate() {
            let range = match od.ref_range(a.start, a.end) {
                Some((x, y)) => format!("{}–{}", group_digits(x as usize), group_digits(y as usize)),
                None => format!("col {}–{}", a.start + 1, a.end),
            };
            let pal = self.palette();
            let swatch = container(Space::new().width(10).height(10)).style(move |_t: &Theme| container::Style {
                background: Some(pal.annotation(a.kind).into()),
                border: iced::Border { radius: 2.0.into(), ..Default::default() },
                ..Default::default()
            });
            anns = anns.push(
                row![
                    swatch,
                    button(column![text(a.name.as_str()).size(12), text(format!("{} · {range}", a.kind)).size(10).style(style::muted)])
                        .padding([2, 4])
                        .style(button::text)
                        .width(Fill)
                        .on_press(Message::GotoAnnotation(i)),
                    button(text("✕").size(11)).padding([2, 6]).style(button::text).on_press(Message::DeleteAnnotation(a.id)),
                ]
                .spacing(6)
                .align_y(Alignment::Center),
            );
        }

        let selection: Element<Message> = match od.selection {
            Some(sel) => {
                let rr = od.ref_range(sel.c0, sel.c1).map(|(a, b)| format!(" · ref {}–{}", group_digits(a as usize), group_digits(b as usize))).unwrap_or_default();
                column![
                    section("Selection"),
                    text(format!("{} columns ({}–{}){rr}", sel.len(), group_digits(sel.c0 + 1), group_digits(sel.c1))).size(12),
                    row![
                        button(text("Oligo report").size(12)).padding([4, 8]).on_press(Message::ShowPopup),
                        tooltip(
                            button(text("Extract slice").size(12))
                                .padding([4, 8])
                                .style(button::secondary)
                                .on_press_maybe((!self.busy()).then_some(Message::ExtractSlice)),
                            container(text("Save these columns of all sequences (without the reference) as a new library document").size(12))
                                .padding(6)
                                .style(container::rounded_box),
                            tooltip::Position::Bottom,
                        ),
                        button(text("Clear").size(12)).padding([4, 8]).style(button::text).on_press(Message::ClearSelection),
                    ]
                    .spacing(6)
                ]
                .spacing(6)
                .into()
            }
            None => column![section("Selection"), text("Drag across the alignment to select columns.").size(11).style(style::muted)].spacing(6).into(),
        };

        let content = column![doc_info, selection, display, highlighting, rows, anns].spacing(18).padding(12);
        container(scrollable(content).height(Fill)).width(OPTIONS_W).height(Fill).style(style::panel).into()
    }

    fn primer_popup<'a>(&'a self, od: &'a OpenDoc, r: &'a PrimerReport) -> Element<'a, Message> {
        let opts = &self.settings.primer;
        let bold = Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT };
        let mono_bold = Font { weight: iced::font::Weight::Bold, ..Font::MONOSPACE };
        let stat = |label: String, n: Option<usize>, pct: Option<f64>| -> Element<'a, Message> {
            row![
                text(label).size(13).width(Length::Fixed(130.0)),
                text(n.map(group_digits).unwrap_or_default()).size(13).width(Length::Fixed(80.0)).align_x(Alignment::End).font(Font::MONOSPACE),
                text(pct.map(|p| format!("{p:.2}%")).unwrap_or_default()).size(13).width(Length::Fixed(80.0)).align_x(Alignment::End).font(Font::MONOSPACE),
            ]
            .into()
        };
        let kv = |label: &'a str, v: String| -> Element<'a, Message> {
            row![text(label).size(13).width(Length::Fixed(130.0)), text(v).size(13).font(Font::MONOSPACE)].into()
        };
        let oligo = ellipsize_middle(&String::from_utf8_lossy(&r.oligo), 60);
        let ref_range = od.ref_range(r.col_start, r.col_end);
        let mut sources = vec![OligoSource::Consensus];
        if od.reference.is_some() {
            sources.push(OligoSource::Reference);
        }
        let source = if od.reference.is_none() { OligoSource::Consensus } else { opts.source };
        let h = &r.mismatch_hist;
        let window = opts.three_prime_window;

        let mut body = column![
            row![
                text("Oligo report").size(15).font(bold),
                Space::new().width(Fill),
                button(text("✕").size(12)).padding([2, 8]).style(button::text).on_press(Message::ClosePopup)
            ]
            .align_y(Alignment::Center),
            text(format!("5'-{oligo}-3'")).size(15).font(mono_bold),
            kv("Length", format!("{}", r.oligo.len())),
            kv(
                "Position",
                match ref_range {
                    Some((a, b)) => format!("ref {}–{} (columns {}–{})", group_digits(a as usize), group_digits(b as usize), r.col_start + 1, r.col_end),
                    None => format!("columns {}–{}", group_digits(r.col_start + 1), group_digits(r.col_end)),
                }
            ),
            kv("Tm", primer::format_tm(r.tm)),
            kv("GC", format!("{:.0}%", r.gc)),
            rule::horizontal(1),
            stat("Sequences".into(), Some(r.sequences), None),
            stat("Perfect match".into(), Some(h[0]), Some(r.percent(h[0]))),
            stat("1 mismatch".into(), Some(h[1]), Some(r.percent(h[1]))),
            stat("2 mismatches".into(), Some(h[2]), Some(r.percent(h[2]))),
            stat("≥3 mismatches".into(), Some(h[3]), Some(r.percent(h[3]))),
            stat(format!("3' mismatch ({window} nt)"), Some(r.three_prime), Some(r.percent(r.three_prime))),
            stat("No coverage".into(), Some(r.no_coverage), None),
        ]
        .spacing(3);
        if r.with_ambiguity > 0 {
            body = body.push(stat("With N / ambiguity".into(), Some(r.with_ambiguity), Some(r.percent(r.with_ambiguity))));
        }
        body = body.push(
            row![
                pick_list(sources, Some(source), Message::PrimerSource).text_size(12).padding(4),
                pick_list(vec![Orientation::Forward, Orientation::Reverse], Some(opts.orientation), Message::PrimerOrientation).text_size(12).padding(4),
            ]
            .spacing(6),
        );
        if !r.per_position.is_empty() && r.sequences > 0 {
            body = body.push(text("Mismatches per position (5'→3')").size(12).style(style::muted));
            body = body.push(
                canvas(ProfileChart { values: &r.per_position, total: r.sequences, window, pal: self.palette(), oligo: &r.oligo })
                    .width(Fill)
                    .height(46),
            );
        }
        if !r.variants.is_empty() {
            body = body.push(text("Most common target sequences").size(12).style(style::muted));
            let mut vars = column![].spacing(1);
            for v in r.variants.iter().take(8) {
                let mm = if v.mismatches == 0 { "✓".to_string() } else { format!("{}mm", v.mismatches) };
                vars = vars.push(row![
                    text(ellipsize_middle(&v.display, 36)).size(12).font(Font::MONOSPACE).width(Fill),
                    text(mm).size(11).font(Font::MONOSPACE).width(Length::Fixed(36.0)).align_x(Alignment::End),
                    text(group_digits(v.count)).size(12).font(Font::MONOSPACE).width(Length::Fixed(64.0)).align_x(Alignment::End),
                    text(format!("{:.1}%", r.percent(v.count))).size(12).font(Font::MONOSPACE).width(Length::Fixed(58.0)).align_x(Alignment::End),
                ]);
            }
            body = body.push(vars);
        }
        if let Some(key) = &r.group_key
            && !r.groups.is_empty() {
                body = body.push(text(format!("Inclusivity by {key}")).size(12).style(style::muted));
                let mut t = column![row![
                    text("Group").size(11).width(Fill).font(bold),
                    text("n").size(11).width(Length::Fixed(52.0)).align_x(Alignment::End).font(bold),
                    text("perfect").size(11).width(Length::Fixed(60.0)).align_x(Alignment::End).font(bold),
                    text("3' mm").size(11).width(Length::Fixed(52.0)).align_x(Alignment::End).font(bold),
                ]]
                .spacing(1);
                for g in r.groups.iter().take(25) {
                    let pct = |n: usize| if g.sequences == 0 { "–".to_string() } else { format!("{:.1}%", 100.0 * n as f64 / g.sequences as f64) };
                    t = t.push(row![
                        text(g.label.as_str()).size(11).width(Fill),
                        text(group_digits(g.sequences)).size(11).width(Length::Fixed(52.0)).align_x(Alignment::End).font(Font::MONOSPACE),
                        text(pct(g.perfect)).size(11).width(Length::Fixed(60.0)).align_x(Alignment::End).font(Font::MONOSPACE),
                        text(pct(g.three_prime)).size(11).width(Length::Fixed(52.0)).align_x(Alignment::End).font(Font::MONOSPACE),
                    ]);
                }
                if r.groups.len() > 25 {
                    t = t.push(text(format!("… {} more groups", r.groups.len() - 25)).size(11).style(style::muted));
                }
                body = body.push(t);
            }
        body = body.push(rule::horizontal(1));
        body = body.push(
            row![
                text_input("Annotation name", &od.annotation_name).on_input(Message::AnnotationName).on_submit(Message::AddAnnotation).padding(4).size(12),
                pick_list(&AnnotationKind::ALL[..], Some(od.annotation_kind), Message::AnnotationKind).text_size(12).padding(4),
            ]
            .spacing(6),
        );
        body = body.push(
            row![
                button(text("Add annotation").size(12)).padding([4, 8]).on_press(Message::AddAnnotation),
                button(text("Extract slice").size(12))
                    .padding([4, 8])
                    .style(button::secondary)
                    .on_press_maybe((!self.busy()).then_some(Message::ExtractSlice)),
                Space::new().width(Fill),
                button(text("Copy oligo").size(12)).padding([4, 8]).style(button::secondary).on_press(Message::CopyOligo),
                button(text("Copy report").size(12)).padding([4, 8]).style(button::secondary).on_press(Message::CopyReport),
            ]
            .spacing(6),
        );
        container(scrollable(body.padding(iced::Padding { right: 10.0, ..Default::default() })).height(Length::Shrink))
            .padding(14)
            .width(Length::Fixed(400.0))
            .max_height(760.0)
            .style(style::card)
            .into()
    }

    fn status_bar(&self) -> Element<'_, Message> {
        let left: Element<Message> = if let Some(e) = &self.error {
            row![
                text(format!("⚠ {e}")).size(12).style(|t: &Theme| text::Style { color: Some(t.extended_palette().danger.base.color) }),
                button(text("Dismiss").size(11)).padding([1, 6]).style(button::text).on_press(Message::Modal(ModalMsg::DismissError)),
            ]
            .spacing(8)
            .align_y(Alignment::Center)
            .into()
        } else {
            text(self.hover_text().unwrap_or_else(|| self.status.clone())).size(12).into()
        };
        let right: Element<Message> = match &self.job {
            Some(j) => {
                let msg = j.progress.message();
                let label = if msg.is_empty() { j.title.clone() } else { msg };
                let mut r = row![text(label).size(12), progress_bar(0.0..=1.0, j.progress.fraction()).length(180).girth(10)]
                    .spacing(8)
                    .align_y(Alignment::Center);
                if j.cancellable {
                    r = r.push(button(text("Cancel").size(11)).padding([1, 8]).style(button::secondary).on_press(Message::CancelJob));
                }
                r.into()
            }
            None => text(self.selection_text().unwrap_or_default()).size(12).into(),
        };
        container(row![left, Space::new().width(Fill), right].align_y(Alignment::Center).spacing(10))
            .padding([4, 12])
            .width(Fill)
            .height(28)
            .style(style::toolbar)
            .into()
    }

    fn hover_text(&self) -> Option<String> {
        let od = self.open.as_ref()?;
        let h = od.hover?;
        let mut s = format!("Column {}", group_digits(h.col + 1));
        if let Some(rc) = &od.ref_coords {
            let p = rc.pos_at.get(h.col).copied().unwrap_or(0);
            s.push_str(&format!(" · Ref {}", group_digits(p as usize)));
            if !rc.is_base[h.col] {
                s.push_str(" (insertion)");
            }
        }
        if let Some(c) = od.consensus.get(h.col).filter(|&&c| c != b' ') {
            s.push_str(&format!(" · Consensus {}", *c as char));
        }
        if let Some(c) = od.conservation.get(h.col).copied().flatten() {
            s.push_str(&format!(" · Conservation {:.2}%", c * 100.0));
        }
        if let Some(id) = od.identity.get(h.col).copied().flatten() {
            s.push_str(&format!(" · Identity {:.1}%", id * 100.0));
        }
        let cov = od.stats.coverage(h.col);
        s.push_str(&format!(" · Coverage {}", group_digits(cov as usize)));
        let row = if h.on_reference { od.reference } else { h.item.and_then(|i| od.item_row(i)) };
        if let Some(r) = row {
            let rw = &od.doc.rows[r];
            let b = rw.at(h.col).map(|b| (b as char).to_string()).unwrap_or_else(|| "no coverage".into());
            s.push_str(&format!(" · {}: {b}", rw.name));
            if let Some(Item::Seq { count, .. }) = h.item.and_then(|i| od.items.get(i))
                && *count > 1 && !h.on_reference {
                    s.push_str(&format!(" (×{count})"));
                }
        }
        Some(s)
    }

    fn selection_text(&self) -> Option<String> {
        let od = self.open.as_ref()?;
        let sel = od.selection?;
        let mut s = format!("Selected {} columns from {} to {}", group_digits(sel.len()), group_digits(sel.c0 + 1), group_digits(sel.c1));
        if let Some((a, b)) = od.ref_range(sel.c0, sel.c1) {
            let n = od.ref_coords.as_ref().map(|rc| (sel.c0..sel.c1).filter(|&c| rc.is_base[c]).count()).unwrap_or(0);
            s.push_str(&format!(" ({} ungapped reference bases from {} to {})", n, group_digits(a as usize), group_digits(b as usize)));
        }
        Some(s)
    }
}

/// Shortens long sequences to `start…end` so they fit on one line.
fn ellipsize_middle(s: &str, max: usize) -> String {
    let n = s.chars().count();
    if n <= max {
        return s.to_string();
    }
    let half = (max - 1) / 2;
    let head: String = s.chars().take(half).collect();
    let tail: String = s.chars().skip(n - half).collect();
    format!("{head}…{tail}")
}

/// Mismatch budget choice for the search pick list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mismatches(u8);

impl std::fmt::Display for Mismatches {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} mm", self.0)
    }
}

/// Bar chart of mismatches per oligo position.
struct ProfileChart<'a> {
    values: &'a [usize],
    total: usize,
    window: usize,
    pal: style::Palette,
    oligo: &'a [u8],
}

impl cv::Program<Message> for ProfileChart<'_> {
    type State = ();

    fn draw(&self, _: &(), renderer: &Renderer, _theme: &Theme, bounds: Rectangle, _cursor: mouse::Cursor) -> Vec<Geometry> {
        let mut f = Frame::new(renderer, bounds.size());
        let n = self.values.len().max(1);
        let w = bounds.width / n as f32;
        let chart_h = bounds.height - 14.0;
        let max = self.values.iter().copied().max().unwrap_or(0).max(1) as f32;
        f.fill_rectangle(Point::new(0.0, chart_h), Size::new(bounds.width, 1.0), self.pal.grid);
        for (i, &v) in self.values.iter().enumerate() {
            let x = i as f32 * w;
            let in_window = i + self.window >= n;
            if in_window {
                f.fill_rectangle(Point::new(x, 0.0), Size::new(w, chart_h), Color { a: 0.08, ..self.pal.id_low });
            }
            if v > 0 {
                // Log-ish scale so rare mismatches remain visible.
                let frac = ((v as f32).ln_1p() / max.ln_1p()).clamp(0.05, 1.0);
                let bh = frac * (chart_h - 2.0);
                let color = if (v as f64) / (self.total as f64) > 0.01 { self.pal.id_low } else { self.pal.id_mid };
                f.fill_rectangle(Point::new(x + 1.0, chart_h - bh), Size::new((w - 2.0).max(1.0), bh), color);
            }
            if w >= 7.0
                && let Some(&b) = self.oligo.get(i) {
                    f.fill_text(cv::Text {
                        content: (b as char).to_string(),
                        position: Point::new(x + w / 2.0, chart_h + 1.0),
                        size: iced::Pixels(10.0),
                        color: self.pal.text_muted,
                        font: Font::MONOSPACE,
                        align_x: iced::widget::text::Alignment::Center,
                        ..cv::Text::default()
                    });
                }
        }
        vec![f.into_geometry()]
    }
}
