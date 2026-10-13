//! Modal dialogs: new/edit sequence, map to reference, export, settings,
//! rename, confirmations.

use super::canvas::group_digits;
use super::state::OpenDoc;
use super::{AppSettings, Message, style};
use crate::fasta::{self, PastedRecord};
use crate::mapper::{MapParams, Sensitivity};
use crate::model::{AnnotationKind, DocInfo, DocKind, RefCoords, Row};
use crate::seq::{self, GAP};
use iced::widget::{Space, button, checkbox, column, container, pick_list, row, rule, text, text_editor, text_input};
use iced::{Element, Fill, Font, Length, Task};

/// Id of the rename text field, focused when the dialog opens.
pub const RENAME_INPUT: &str = "rename-input";

/// Focuses the rename field and selects its text.
pub fn focus_rename<T: Send + 'static>() -> Task<T> {
    Task::batch([iced::widget::operation::focus(RENAME_INPUT), iced::widget::operation::select_all(RENAME_INPUT)])
}

/// Id of the sequence text box of the new/edit sequence dialog.
pub const SEQUENCE_INPUT: &str = "sequence-input";

/// Focuses the sequence text box (ready for pasting).
pub fn focus_sequence<T: Send + 'static>() -> Task<T> {
    iced::widget::operation::focus(SEQUENCE_INPUT)
}

/// Id of the annotation name field, focused when the dialog opens.
pub const ANNOTATION_NAME: &str = "annotation-name";

/// Focuses the annotation name and selects its text.
pub fn focus_annotation<T: Send + 'static>() -> Task<T> {
    Task::batch([iced::widget::operation::focus(ANNOTATION_NAME), iced::widget::operation::select_all(ANNOTATION_NAME)])
}

pub enum Modal {
    Sequence(SequenceDialog),
    Annotation(AnnotationDialog),
    ConfirmDeleteAnnotation { id: i64, name: String, detail: String },
    Map(MapDialog),
    Export(ExportDialog),
    Settings(SettingsDialog),
    Rename { id: i64, text: String, folder: bool },
    ConfirmDelete { id: i64, name: String, detail: String },
    ConfirmDeleteFolder { id: i64, name: String, documents: usize, subfolders: usize },
    NameFields { delimiter: String, keys: String },
    Info { title: String, text: String },
}

#[derive(Debug, Clone)]
pub enum ModalMsg {
    RenameInput(String),
    RenameConfirm,
    DelimiterInput(String),
    KeysInput(String),
    ApplyNameFields,
    DismissError,
}

// ---------------------------------------------------------------- new / edit sequence

/// The single-sequence document edited by the sequence dialog.
pub struct EditTarget {
    pub id: i64,
    pub doc_name: String,
    pub kind: DocKind,
    /// The stored sequence (name, description, metadata, placement).
    pub row: Row,
    /// End columns of the document's annotations.
    pub annotation_ends: Vec<usize>,
}

impl EditTarget {
    /// The stored sequence with its leading gaps, as shown in the text box.
    fn gapped(&self) -> Vec<u8> {
        let mut v = vec![GAP; self.row.start];
        v.extend_from_slice(&self.row.data);
        v
    }
}

/// Creates a document from typed or pasted text (plain sequence or FASTA), or
/// edits the sequence of a single-sequence document.
pub struct SequenceDialog {
    /// `None`: create a new document.
    pub edit: Option<EditTarget>,
    /// Folder that receives a new document.
    pub folder: Option<i64>,
    pub name: String,
    pub description: String,
    pub content: text_editor::Content,
    /// The parsed text, refreshed on every edit.
    pub records: Vec<PastedRecord>,
    /// Name and description last suggested (FASTA header or stored values);
    /// they follow the text until the user edits the fields.
    suggested: (String, String),
    /// Summary of the parsed text and notes, computed on edits (not per frame).
    summary: String,
    notes: Vec<String>,
    /// Record names and lengths when several records were entered.
    listing: Vec<String>,
}

#[derive(Debug, Clone)]
pub enum SequenceMsg {
    Name(String),
    Description(String),
    Edit(text_editor::Action),
}

/// Bases per line when an existing sequence is shown for editing.
const EDIT_LINE: usize = 60;

impl SequenceDialog {
    pub fn create(folder: Option<i64>) -> SequenceDialog {
        Self::with(None, folder, String::new())
    }

    pub fn edit(target: EditTarget) -> SequenceDialog {
        let text = target.gapped().chunks(EDIT_LINE).map(String::from_utf8_lossy).collect::<Vec<_>>().join("\n");
        Self::with(Some(target), None, text)
    }

    fn with(edit: Option<EditTarget>, folder: Option<i64>, text: String) -> SequenceDialog {
        let mut d = SequenceDialog {
            edit,
            folder,
            name: String::new(),
            description: String::new(),
            content: text_editor::Content::with_text(&text),
            records: Vec::new(),
            suggested: (String::new(), String::new()),
            summary: String::new(),
            notes: Vec::new(),
            listing: Vec::new(),
        };
        d.suggested = d.base_suggestion();
        (d.name, d.description) = d.suggested.clone();
        d.refresh();
        d
    }

    /// Name and description without a FASTA header: empty, or the stored ones.
    fn base_suggestion(&self) -> (String, String) {
        match &self.edit {
            Some(t) => (t.doc_name.clone(), t.row.description.clone()),
            None => (String::new(), String::new()),
        }
    }

    pub fn update(&mut self, m: SequenceMsg) {
        match m {
            SequenceMsg::Name(s) => self.name = s,
            SequenceMsg::Description(s) => self.description = s,
            SequenceMsg::Edit(action) => {
                let edit = action.is_edit();
                self.content.perform(action);
                if edit {
                    self.refresh();
                }
            }
        }
    }

    /// Re-parses the text; a single FASTA header fills in name and description
    /// unless the user changed them.
    fn refresh(&mut self) {
        self.records = fasta::parse_pasted(&self.content.text());
        let suggestion = match self.records.as_slice() {
            [r] if r.id.is_some() => (r.id.clone().unwrap_or_default(), r.description.clone()),
            _ => self.base_suggestion(),
        };
        if self.name == self.suggested.0 {
            self.name = suggestion.0.clone();
        }
        if self.description == self.suggested.1 {
            self.description = suggestion.1.clone();
        }
        self.suggested = suggestion;

        self.notes.clear();
        self.listing.clear();
        let unknown: usize = self.records.iter().map(|r| r.unknown).sum();
        self.summary = match self.records.as_slice() {
            [] => String::new(),
            [r] => {
                let bases = r.seq.iter().filter(|&&b| b != GAP).count();
                let gaps = r.seq.len() - bases;
                let ambiguous = r.seq.iter().filter(|&&b| b != GAP && !seq::is_acgt(b)).count();
                let mut s = format!("{} bp · GC {:.1}%", group_digits(bases), seq::gc_percent(&r.seq));
                if ambiguous > 0 {
                    s.push_str(&format!(" · {} ambiguous (N, R, Y…)", group_digits(ambiguous)));
                }
                if gaps > 0 {
                    s.push_str(&format!(" · {} gaps", group_digits(gaps)));
                }
                s
            }
            recs => {
                let lens: Vec<usize> = recs.iter().map(|r| r.seq.iter().filter(|&&b| b != GAP).count()).collect();
                let (min, max) = (lens.iter().min().copied().unwrap_or(0), lens.iter().max().copied().unwrap_or(0));
                let range = if min == max { format!("{} bp each", group_digits(min)) } else { format!("{}–{} bp", group_digits(min), group_digits(max)) };
                const SHOWN: usize = 6;
                for (i, (r, len)) in recs.iter().zip(&lens).take(SHOWN).enumerate() {
                    let name = r.id.clone().unwrap_or_else(|| format!("Sequence {}", i + 1));
                    self.listing.push(format!("{name} · {} bp", group_digits(*len)));
                }
                if recs.len() > SHOWN {
                    self.listing.push(format!("… {} more", recs.len() - SHOWN));
                }
                format!("{} sequences · {range} → new {}", recs.len(), fasta::pasted_kind(recs).label().to_lowercase())
            }
        };
        if unknown > 0 {
            self.notes.push(format!("{} character(s) that are not nucleotide codes will be stored as N.", group_digits(unknown)));
        }
        if let (Some(t), [r]) = (&self.edit, self.records.as_slice())
            && !t.annotation_ends.is_empty()
            && r.seq != t.gapped()
        {
            let end = r.seq.iter().rposition(|&b| b != GAP).map_or(0, |p| p + 1);
            let n = t.annotation_ends.len();
            let dropped = t.annotation_ends.iter().filter(|&&e| e > end).count();
            let mut note = if n == 1 {
                "This document has 1 annotation; its position is not shifted when bases are inserted or removed.".to_string()
            } else {
                format!("This document has {n} annotations; their positions are not shifted when bases are inserted or removed.")
            };
            note.push_str(&match (n, dropped) {
                (_, 0) => String::new(),
                (1, _) => " It extends past the new end and will be removed.".into(),
                (_, 1) => " 1 of them extends past the new end and will be removed.".into(),
                (_, d) => format!(" {d} of them extend past the new end and will be removed."),
            });
            self.notes.push(note);
        }
    }

    /// Why the dialog cannot be confirmed yet.
    pub fn problem(&self) -> Option<String> {
        if self.records.is_empty() {
            return Some("Type or paste a sequence.".into());
        }
        if let Some((i, r)) = self.records.iter().enumerate().find(|(_, r)| r.seq.iter().all(|&b| b == GAP)) {
            let name = r.id.clone().unwrap_or_else(|| format!("Sequence {}", i + 1));
            return Some(if self.records.len() == 1 { "The record has no sequence yet.".into() } else { format!("“{name}” has no sequence.") });
        }
        if self.edit.is_some() && self.records.len() > 1 {
            return Some(format!(
                "{} records found: this document holds one sequence. Use File → New sequence… to create a sequence list.",
                self.records.len()
            ));
        }
        if self.name.trim().is_empty() {
            return Some("Enter a name.".into());
        }
        None
    }

    pub fn view(&self) -> Element<'_, Message> {
        let s = Message::Sequence;
        let multi = self.records.len() > 1;
        let label = |t: &'static str| text(t).size(13).width(Length::Fixed(96.0));
        let name_hint = match (&self.edit, multi) {
            (None, true) => "Name of the new document",
            _ => "Name (or taken from a FASTA header)",
        };
        // The widget structure stays the same in every mode: a changed layout
        // in front of the text box would reset it and drop its focus.
        let list = multi && self.edit.is_none();
        let description = if list {
            text_input("Taken from each FASTA header", "").padding(6)
        } else {
            text_input("Optional; key=value pairs become metadata", &self.description)
                .on_input(move |v| s(SequenceMsg::Description(v)))
                .on_submit(Message::SaveSequence)
                .padding(6)
        };
        let name = text_input(name_hint, &self.name).on_input(move |v| s(SequenceMsg::Name(v))).on_submit(Message::SaveSequence).padding(6);
        let mut body = column![
            row![label(if list { "Document" } else { "Name" }), name].spacing(8).align_y(iced::Center),
            row![label("Description"), description].spacing(8).align_y(iced::Center),
        ]
        .spacing(10);
        body = body.push(
            text_editor(&self.content)
                .id(SEQUENCE_INPUT)
                .placeholder("Paste or type a nucleotide sequence, or FASTA text with one or more records")
                .on_action(move |a| s(SequenceMsg::Edit(a)))
                .font(Font::MONOSPACE)
                .size(13)
                .padding(8)
                .height(Length::Fixed(300.0))
                .wrapping(text::Wrapping::WordOrGlyph),
        );
        let help = if self.edit.is_some() {
            "Spaces, line breaks and numbers are ignored; U becomes T; other characters that are not nucleotide codes become N."
        } else {
            "Spaces, line breaks and numbers (GenBank listings) are ignored; U becomes T; other characters that are not nucleotide codes become N. \
             A FASTA header (>name description) fills in the name; several records create a sequence list."
        };
        body = body.push(text(help).size(12).style(style::muted));
        // Always one summary line, so the dialog does not jump on the first key.
        body = body.push(text(if self.summary.is_empty() { "0 bp" } else { self.summary.as_str() }).size(13).font(Font {
            weight: iced::font::Weight::Semibold,
            ..Font::DEFAULT
        }));
        if !self.listing.is_empty() {
            body = body.push(column(self.listing.iter().map(|l| text(l.as_str()).size(12).font(Font::MONOSPACE).into())).spacing(1));
        }
        for n in &self.notes {
            body = body.push(text(n.as_str()).size(12).style(|t: &iced::Theme| text::Style { color: Some(t.extended_palette().danger.base.color) }));
        }
        let problem = self.problem();
        let confirm = if self.edit.is_some() { "Save" } else { "Create" };
        let buttons = row![
            text(problem.clone().unwrap_or_default()).size(12).style(style::muted).width(Fill),
            button("Cancel").on_press(Message::CloseModal).style(button::secondary),
            button(confirm).on_press_maybe(problem.is_none().then_some(Message::SaveSequence)).style(button::primary),
        ]
        .spacing(10)
        .align_y(iced::Center);
        let title = if self.edit.is_some() { "Edit sequence" } else { "New sequence" };
        dialog_frame(title, body.into(), buttons.into(), 660.0)
    }
}

// ---------------------------------------------------------------- annotation properties

/// Creates an annotation from the selection or edits one: name, type, range.
pub struct AnnotationDialog {
    /// `None`: a new annotation.
    pub id: Option<i64>,
    pub name: String,
    pub kind: AnnotationKind,
    pub start: String,
    pub end: String,
    /// Columns and field texts when the dialog opened. A field the user does
    /// not change keeps its column (reference positions cannot address
    /// insertion columns).
    initial: (usize, usize, String, String),
    /// Reference coordinates of the range fields; `None`: alignment columns.
    coords: Option<RefCoords>,
    width: usize,
}

#[derive(Debug, Clone)]
pub enum AnnotationMsg {
    Name(String),
    Kind(AnnotationKind),
    Start(String),
    End(String),
}

/// Range fields for columns `[c0, c1)`: 1-based reference positions, or columns.
fn range_text(c0: usize, c1: usize, coords: Option<&RefCoords>) -> (String, String) {
    match coords.and_then(|rc| rc.range(c0, c1)) {
        Some((a, b)) => (a.to_string(), b.to_string()),
        None => ((c0 + 1).to_string(), c1.to_string()),
    }
}

impl AnnotationDialog {
    pub fn new(id: Option<i64>, name: String, kind: AnnotationKind, c0: usize, c1: usize, od: &OpenDoc) -> AnnotationDialog {
        // Reference positions when the range holds reference bases.
        let coords = od.ref_coords.clone().filter(|rc| rc.range(c0, c1).is_some());
        let (start, end) = range_text(c0, c1, coords.as_ref());
        AnnotationDialog {
            id,
            name,
            kind,
            initial: (c0, c1, start.clone(), end.clone()),
            start,
            end,
            coords,
            width: od.width(),
        }
    }

    pub fn update(&mut self, m: AnnotationMsg) {
        match m {
            AnnotationMsg::Name(s) => self.name = s,
            AnnotationMsg::Kind(k) => self.kind = k,
            AnnotationMsg::Start(s) => self.start = s,
            AnnotationMsg::End(s) => self.end = s,
        }
    }

    /// The range as columns `[c0, c1)`, or why the fields are invalid.
    pub fn range(&self) -> Result<(usize, usize), String> {
        let parse = |s: &str, what: &str| s.trim().replace(',', "").parse::<usize>().map_err(|_| format!("{what} must be a whole number."));
        let (p0, p1) = (parse(&self.start, "Start")?, parse(&self.end, "End")?);
        let max = self.coords.as_ref().map_or(self.width, |rc| rc.col_of.len());
        if p0 < 1 || p1 > max {
            return Err(format!("The range must lie within 1–{}.", group_digits(max)));
        }
        if p1 < p0 {
            return Err("The end must not be before the start.".into());
        }
        let (c0, c1) = match &self.coords {
            Some(rc) => rc.columns(p0, p1).ok_or("Invalid range.")?,
            None => (p0 - 1, p1),
        };
        let (i0, i1, s0, s1) = &self.initial;
        let c0 = if self.start.trim() == s0 { *i0 } else { c0 };
        let c1 = if self.end.trim() == s1 { *i1 } else { c1 };
        if c1 <= c0 {
            return Err("The end must not be before the start.".into());
        }
        Ok((c0, c1))
    }

    pub fn problem(&self) -> Option<String> {
        if self.name.trim().is_empty() {
            return Some("Enter a name.".into());
        }
        self.range().err()
    }

    pub fn view(&self) -> Element<'_, Message> {
        let a = Message::AnnotationEdit;
        let label = |t: &'static str| text(t).size(13).width(Length::Fixed(96.0));
        let unit = if self.coords.is_some() { "reference position" } else { "alignment column" };
        let field = |value: &str, f: fn(String) -> AnnotationMsg| {
            text_input("", value).on_input(move |s| a(f(s))).on_submit(Message::SaveAnnotation).padding(6).width(110)
        };
        let info: Element<Message> = match self.range() {
            Ok((c0, c1)) => {
                let cols = format!("{} column(s), {}–{}", group_digits(c1 - c0), group_digits(c0 + 1), group_digits(c1));
                let s = match &self.coords {
                    Some(rc) => {
                        let bases = (c0..c1).filter(|&c| rc.is_base[c]).count();
                        format!("{} reference bases · {cols}", group_digits(bases))
                    }
                    None => cols,
                };
                text(s).size(12).style(style::muted).into()
            }
            Err(e) => text(e).size(12).style(|t: &iced::Theme| text::Style { color: Some(t.extended_palette().danger.base.color) }).into(),
        };
        let body = column![
            row![
                label("Name"),
                text_input("Annotation name", &self.name)
                    .id(ANNOTATION_NAME)
                    .on_input(move |s| a(AnnotationMsg::Name(s)))
                    .on_submit(Message::SaveAnnotation)
                    .padding(6)
            ]
            .spacing(8)
            .align_y(iced::Center),
            row![label("Type"), pick_list(&AnnotationKind::ALL[..], Some(self.kind), move |k| a(AnnotationMsg::Kind(k))).width(Fill)]
                .spacing(8)
                .align_y(iced::Center),
            row![label("Range"), field(&self.start, AnnotationMsg::Start), text("to").size(13), field(&self.end, AnnotationMsg::End), text(unit).size(12).style(style::muted)]
                .spacing(8)
                .align_y(iced::Center),
            row![Space::new().width(Length::Fixed(104.0)), info],
        ]
        .spacing(10);
        let mut buttons = row![].spacing(10).align_y(iced::Center);
        if let Some(id) = self.id {
            buttons = buttons.push(button("Delete…").on_press(Message::AskDeleteAnnotation(id)).style(button::danger));
        }
        let ok = self.problem().is_none();
        buttons = buttons.push(Space::new().width(Fill)).push(button("Cancel").on_press(Message::CloseModal).style(button::secondary)).push(
            button(if self.id.is_some() { "Save" } else { "Add" }).on_press_maybe(ok.then_some(Message::SaveAnnotation)).style(button::primary),
        );
        dialog_frame(if self.id.is_some() { "Annotation properties" } else { "New annotation" }, body.into(), buttons.into(), 520.0)
    }
}

// ---------------------------------------------------------------- mapping

#[derive(Debug, Clone, PartialEq)]
pub enum RefChoice {
    ThisDocument(Option<usize>),
    OtherDocument(i64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct RefOption {
    pub row: usize,
    pub name: String,
}

impl std::fmt::Display for RefOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.name)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DocOption {
    pub id: i64,
    pub name: String,
}

impl std::fmt::Display for DocOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.name)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FineTuning(pub usize);

impl std::fmt::Display for FineTuning {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            0 => f.write_str("None (single pass)"),
            n => write!(f, "Iterate up to {n} time{}", if n > 1 { "s" } else { "" }),
        }
    }
}

pub struct MapDialog {
    pub reference: RefChoice,
    pub ref_filter: String,
    pub candidates: Vec<RefOption>,
    pub docs: Vec<DocOption>,
    pub sensitivity: Sensitivity,
    pub word_len: String,
    pub repeats: String,
    pub min_identity: String,
    pub min_aligned: String,
    pub fine_tuning: usize,
    pub both_strands: bool,
    pub trim: bool,
    pub name: String,
    pub only_selected: bool,
    pub save_unmapped: bool,
    doc_name: String,
    n_rows: usize,
    n_selected: usize,
}

#[derive(Debug, Clone)]
pub enum MapMsg {
    UseThisDocument(bool),
    RefFilter(String),
    RefRow(RefOption),
    RefDoc(DocOption),
    Sensitivity(Sensitivity),
    WordLen(String),
    Repeats(String),
    MinIdentity(String),
    MinAligned(String),
    FineTuning(FineTuning),
    BothStrands(bool),
    Trim(bool),
    Name(String),
    OnlySelected(bool),
    SaveUnmapped(bool),
}

const MAX_CANDIDATES: usize = 300;

impl MapDialog {
    pub fn new(od: &OpenDoc, last: &MapParams, library: &[DocInfo]) -> MapDialog {
        let selected_single = (od.selected_rows.len() == 1).then(|| *od.selected_rows.iter().next().unwrap());
        let reference = od.reference.or(selected_single).or((!od.doc.rows.is_empty()).then_some(0));
        let mut d = MapDialog {
            reference: RefChoice::ThisDocument(reference),
            ref_filter: String::new(),
            candidates: Vec::new(),
            docs: library.iter().filter(|x| x.id != od.info.id).map(|x| DocOption { id: x.id, name: x.name.clone() }).collect(),
            sensitivity: last.sensitivity,
            word_len: last.word_len.to_string(),
            repeats: last.max_word_repeats.to_string(),
            min_identity: format!("{}", last.min_identity),
            min_aligned: last.min_aligned.to_string(),
            fine_tuning: last.fine_tuning,
            both_strands: last.both_strands,
            trim: last.trim_to_reference,
            name: String::new(),
            only_selected: false,
            save_unmapped: true,
            doc_name: od.info.name.clone(),
            n_rows: od.doc.rows.len(),
            n_selected: od.selected_rows.len(),
        };
        d.refresh_candidates(od);
        d.name = d.default_name(od);
        d
    }

    fn default_name(&self, od: &OpenDoc) -> String {
        let ref_name = match &self.reference {
            RefChoice::ThisDocument(Some(r)) => od.doc.rows[*r].name.clone(),
            RefChoice::OtherDocument(id) => self.docs.iter().find(|d| d.id == *id).map(|d| d.name.clone()).unwrap_or_default(),
            RefChoice::ThisDocument(None) => String::new(),
        };
        format!("{} mapped to {}", self.doc_name, ref_name)
    }

    fn refresh_candidates(&mut self, od: &OpenDoc) {
        let f = self.ref_filter.to_lowercase();
        self.candidates = od
            .doc
            .rows
            .iter()
            .enumerate()
            .filter(|(_, r)| f.is_empty() || r.name.to_lowercase().contains(&f))
            .take(MAX_CANDIDATES)
            .map(|(i, r)| RefOption { row: i, name: r.name.clone() })
            .collect();
    }

    pub fn update(&mut self, m: MapMsg, od: Option<&OpenDoc>) {
        let old_default = od.map(|o| self.default_name(o));
        match m {
            MapMsg::UseThisDocument(b) => {
                self.reference = if b {
                    RefChoice::ThisDocument(od.and_then(|o| o.reference).or(Some(0)))
                } else {
                    match self.docs.first() {
                        Some(d) => RefChoice::OtherDocument(d.id),
                        None => self.reference.clone(),
                    }
                }
            }
            MapMsg::RefFilter(s) => {
                self.ref_filter = s;
                if let Some(od) = od {
                    self.refresh_candidates(od);
                    if self.candidates.len() == 1 {
                        self.reference = RefChoice::ThisDocument(Some(self.candidates[0].row));
                    }
                }
            }
            MapMsg::RefRow(r) => self.reference = RefChoice::ThisDocument(Some(r.row)),
            MapMsg::RefDoc(d) => self.reference = RefChoice::OtherDocument(d.id),
            MapMsg::Sensitivity(s) => {
                self.sensitivity = s;
                if let Some((k, rep, id)) = s.preset() {
                    self.word_len = k.to_string();
                    self.repeats = rep.to_string();
                    self.min_identity = format!("{id}");
                }
            }
            MapMsg::WordLen(s) => {
                self.word_len = s;
                self.sensitivity = Sensitivity::Custom;
            }
            MapMsg::Repeats(s) => {
                self.repeats = s;
                self.sensitivity = Sensitivity::Custom;
            }
            MapMsg::MinIdentity(s) => {
                self.min_identity = s;
                self.sensitivity = Sensitivity::Custom;
            }
            MapMsg::MinAligned(s) => self.min_aligned = s,
            MapMsg::FineTuning(f) => self.fine_tuning = f.0,
            MapMsg::BothStrands(b) => self.both_strands = b,
            MapMsg::Trim(b) => self.trim = b,
            MapMsg::Name(s) => self.name = s,
            MapMsg::OnlySelected(b) => self.only_selected = b,
            MapMsg::SaveUnmapped(b) => self.save_unmapped = b,
        }
        // Keep the suggested name in sync unless the user edited it.
        if let (Some(od), Some(old)) = (od, old_default)
            && self.name == old {
                self.name = self.default_name(od);
            }
    }

    pub fn params(&self) -> Result<MapParams, String> {
        let int = |s: &str, what: &str| s.trim().parse::<usize>().map_err(|_| format!("{what} must be a whole number"));
        let word_len = int(&self.word_len, "Word length")?;
        if !(6..=31).contains(&word_len) {
            return Err("Word length must be between 6 and 31".into());
        }
        let min_identity: f32 = self.min_identity.trim().trim_end_matches('%').parse().map_err(|_| "Minimum identity must be a number".to_string())?;
        Ok(MapParams {
            sensitivity: self.sensitivity,
            word_len,
            max_word_repeats: int(&self.repeats, "Repeat limit")?,
            min_identity: min_identity.clamp(0.0, 100.0),
            min_aligned: int(&self.min_aligned, "Minimum aligned length")?,
            fine_tuning: self.fine_tuning,
            both_strands: self.both_strands,
            trim_to_reference: self.trim,
        })
    }

    pub fn view(&self) -> Element<'_, Message> {
        let m = |x: MapMsg| Message::Map(x);
        let this_doc = matches!(self.reference, RefChoice::ThisDocument(_));
        let ref_section: Element<Message> = match &self.reference {
            RefChoice::ThisDocument(sel) => {
                let selected = sel.and_then(|r| self.candidates.iter().find(|c| c.row == r).cloned());
                column![
                    text_input("Filter sequence names…", &self.ref_filter).on_input(move |s| m(MapMsg::RefFilter(s))).padding(6),
                    pick_list(self.candidates.clone(), selected, move |r| m(MapMsg::RefRow(r)))
                        .placeholder("Choose the reference sequence")
                        .width(Fill),
                ]
                .spacing(6)
                .into()
            }
            RefChoice::OtherDocument(id) => {
                let selected = self.docs.iter().find(|d| d.id == *id).cloned();
                column![
                    pick_list(self.docs.clone(), selected, move |d| m(MapMsg::RefDoc(d))).width(Fill),
                    text("Uses that document's reference sequence (or its first sequence).").size(12).style(style::muted),
                ]
                .spacing(6)
                .into()
            }
        };
        let custom = self.sensitivity == Sensitivity::Custom;
        let field = |label: &'static str, value: &str, f: fn(String) -> MapMsg| {
            row![text(label).width(Fill), text_input("", value).on_input(move |s| Message::Map(f(s))).padding(5).width(80)]
                .spacing(8)
                .align_y(iced::Center)
        };
        let queries = match (&self.reference, self.only_selected && self.n_selected > 0) {
            (_, true) => self.n_selected,
            (RefChoice::OtherDocument(_), false) => self.n_rows,
            (RefChoice::ThisDocument(_), false) => self.n_rows.saturating_sub(1),
        };
        let left = column![
            text("Sequences").size(15),
            text(format!("Maps the sequences of “{}” (about {} sequences).", self.doc_name, super::canvas::group_digits(queries))).size(13),
            checkbox(self.only_selected)
                .label(format!("Only the {} selected rows", self.n_selected))
                .on_toggle_maybe((self.n_selected > 0).then_some(move |b| m(MapMsg::OnlySelected(b)))),
            rule::horizontal(1),
            text("Reference").size(15),
            checkbox(this_doc).label("A sequence of this document").on_toggle(move |b| m(MapMsg::UseThisDocument(b))),
            checkbox(!this_doc)
                .label("Another document from the library")
                .on_toggle_maybe((!self.docs.is_empty()).then_some(move |b: bool| m(MapMsg::UseThisDocument(!b)))),
            ref_section,
            rule::horizontal(1),
            text("Result").size(15),
            text_input("Name of the contig", &self.name).on_input(move |s| m(MapMsg::Name(s))).padding(6),
            checkbox(self.save_unmapped).label("Save unmapped sequences as a separate document").on_toggle(move |b| m(MapMsg::SaveUnmapped(b))),
        ]
        .spacing(10)
        .width(Fill);
        let right = column![
            text("Sensitivity").size(15),
            pick_list(Sensitivity::ALL, Some(self.sensitivity), move |s| m(MapMsg::Sensitivity(s))).width(Fill),
            field("Index word length (k)", &self.word_len, MapMsg::WordLen),
            field("Skip words repeated more than", &self.repeats, MapMsg::Repeats),
            field("Minimum identity (%)", &self.min_identity, MapMsg::MinIdentity),
            field("Minimum aligned length (bp)", &self.min_aligned, MapMsg::MinAligned),
            if custom { text("Custom values.").size(12).style(style::muted) } else { text("Editing a value switches to custom sensitivity.").size(12).style(style::muted) },
            rule::horizontal(1),
            text("Fine tuning").size(15),
            pick_list((0..=5).map(FineTuning).collect::<Vec<_>>(), Some(FineTuning(self.fine_tuning)), move |f| m(MapMsg::FineTuning(f))).width(Fill),
            text("Re-maps all sequences to the consensus of the previous round so indels line up between sequences.")
                .size(12)
                .style(style::muted),
            checkbox(self.both_strands).label("Map both orientations (auto reverse complement)").on_toggle(move |b| m(MapMsg::BothStrands(b))),
            checkbox(self.trim).label("Trim to the reference (clip overhanging ends)").on_toggle(move |b| m(MapMsg::Trim(b))),
        ]
        .spacing(10)
        .width(Fill);
        let body = row![left, right].spacing(32);
        let buttons = row![
            Space::new().width(Fill),
            button("Cancel").on_press(Message::CloseModal).style(button::secondary),
            button("Map to reference").on_press(Message::RunMapping).style(button::primary),
        ]
        .spacing(10);
        dialog_frame("Map to Reference", body.into(), buttons.into(), 900.0)
    }
}

// ---------------------------------------------------------------- export

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportRows {
    All,
    Shown,
    Selected,
}

impl std::fmt::Display for ExportRows {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            ExportRows::All => "All sequences",
            ExportRows::Shown => "Sequences shown (filter / folded groups applied)",
            ExportRows::Selected => "Selected rows",
        })
    }
}

pub struct ExportDialog {
    pub rows: ExportRows,
    pub selection_only: bool,
    pub gapped: bool,
    pub include_reference: bool,
    pub include_consensus: bool,
    /// Leave out the columns where the reference has a gap.
    pub omit_ref_gaps: bool,
    has_selection: bool,
    has_reference: bool,
    selection_len: usize,
    n_selected: usize,
    doc_name: String,
}

#[derive(Debug, Clone)]
pub enum ExportMsg {
    Rows(ExportRows),
    SelectionOnly(bool),
    Gapped(bool),
    IncludeReference(bool),
    IncludeConsensus(bool),
    OmitRefGaps(bool),
    Save,
}

impl ExportDialog {
    pub fn new(od: &OpenDoc) -> ExportDialog {
        ExportDialog {
            rows: if od.selected_rows.is_empty() { ExportRows::All } else { ExportRows::Selected },
            selection_only: false,
            gapped: od.info.kind != crate::model::DocKind::Sequences,
            include_reference: od.reference.is_some(),
            include_consensus: false,
            omit_ref_gaps: od.cols.is_subset(),
            has_selection: od.selection.is_some(),
            has_reference: od.reference.is_some(),
            selection_len: od.selection.map(|s| s.len()).unwrap_or(0),
            n_selected: od.selected_rows.len(),
            doc_name: od.info.name.clone(),
        }
    }

    pub fn update(&mut self, m: ExportMsg) -> Option<Task<Message>> {
        match m {
            ExportMsg::Rows(r) => self.rows = r,
            ExportMsg::SelectionOnly(b) => self.selection_only = b,
            ExportMsg::Gapped(b) => self.gapped = b,
            ExportMsg::IncludeReference(b) => self.include_reference = b,
            ExportMsg::IncludeConsensus(b) => self.include_consensus = b,
            ExportMsg::OmitRefGaps(b) => self.omit_ref_gaps = b,
            ExportMsg::Save => {
                let file_name = format!("{}.fasta", sanitize(&self.doc_name));
                return Some(Task::perform(
                    rfd::AsyncFileDialog::new()
                        .set_title("Export FASTA")
                        .set_file_name(file_name)
                        .add_filter("FASTA", &["fasta", "fa", "fas"])
                        .save_file(),
                    |f| Message::ExportPicked(f.map(|f| f.path().to_path_buf())),
                ));
            }
        }
        None
    }

    pub fn view(&self) -> Element<'_, Message> {
        let e = Message::Export;
        let mut options = vec![ExportRows::All, ExportRows::Shown];
        if self.n_selected > 0 {
            options.push(ExportRows::Selected);
        }
        let body = column![
            text("Sequences").size(15),
            pick_list(options, Some(self.rows), move |r| e(ExportMsg::Rows(r))).width(Fill),
            checkbox(self.selection_only)
                .label(format!("Only the selected columns ({} columns)", self.selection_len))
                .on_toggle_maybe(self.has_selection.then_some(move |b| e(ExportMsg::SelectionOnly(b)))),
            checkbox(self.gapped).label("Keep alignment gaps (aligned FASTA)").on_toggle(move |b| e(ExportMsg::Gapped(b))),
            checkbox(self.include_reference)
                .label("Include the reference sequence")
                .on_toggle_maybe(self.has_reference.then_some(move |b| e(ExportMsg::IncludeReference(b)))),
            checkbox(self.include_consensus).label("Include the consensus sequence").on_toggle(move |b| e(ExportMsg::IncludeConsensus(b))),
            checkbox(self.omit_ref_gaps)
                .label("Leave out columns with gaps in the reference")
                .on_toggle_maybe(self.has_reference.then_some(move |b| e(ExportMsg::OmitRefGaps(b)))),
        ]
        .spacing(10);
        let buttons = row![
            Space::new().width(Fill),
            button("Cancel").on_press(Message::CloseModal).style(button::secondary),
            button("Save as…").on_press(e(ExportMsg::Save)).style(button::primary),
        ]
        .spacing(10);
        dialog_frame("Export FASTA", body.into(), buttons.into(), 480.0)
    }
}

fn sanitize(s: &str) -> String {
    s.chars().map(|c| if c.is_alphanumeric() || " -_.()".contains(c) { c } else { '_' }).collect()
}

// ---------------------------------------------------------------- settings

pub struct SettingsDialog {
    oligo: String,
    na: String,
    mg: String,
    dntp: String,
    window: String,
    dark: bool,
}

#[derive(Debug, Clone)]
pub enum SettingsMsg {
    Oligo(String),
    Na(String),
    Mg(String),
    Dntp(String),
    Window(String),
    Dark(bool),
    Defaults,
}

impl SettingsDialog {
    pub fn new(s: &AppSettings) -> SettingsDialog {
        let tm = s.primer.tm;
        SettingsDialog {
            oligo: format!("{}", tm.oligo_nm),
            na: format!("{}", tm.monovalent_mm),
            mg: format!("{}", tm.divalent_mm),
            dntp: format!("{}", tm.dntp_mm),
            window: s.primer.three_prime_window.to_string(),
            dark: s.dark,
        }
    }

    pub fn update(&mut self, m: SettingsMsg) {
        match m {
            SettingsMsg::Oligo(s) => self.oligo = s,
            SettingsMsg::Na(s) => self.na = s,
            SettingsMsg::Mg(s) => self.mg = s,
            SettingsMsg::Dntp(s) => self.dntp = s,
            SettingsMsg::Window(s) => self.window = s,
            SettingsMsg::Dark(b) => self.dark = b,
            SettingsMsg::Defaults => {
                let d = AppSettings { dark: self.dark, ..AppSettings::default() };
                *self = SettingsDialog::new(&d);
            }
        }
    }

    pub fn apply(&self, s: &mut AppSettings) {
        let num = |v: &str, old: f64| v.trim().replace(',', ".").parse::<f64>().ok().filter(|x| *x >= 0.0).unwrap_or(old);
        let tm = &mut s.primer.tm;
        tm.oligo_nm = num(&self.oligo, tm.oligo_nm).max(0.001);
        tm.monovalent_mm = num(&self.na, tm.monovalent_mm);
        tm.divalent_mm = num(&self.mg, tm.divalent_mm);
        tm.dntp_mm = num(&self.dntp, tm.dntp_mm);
        if let Ok(w) = self.window.trim().parse::<usize>() {
            s.primer.three_prime_window = w.clamp(1, 20);
        }
        s.dark = self.dark;
    }

    pub fn view(&self) -> Element<'_, Message> {
        let field = |label: &'static str, value: &str, f: fn(String) -> SettingsMsg| {
            row![text(label).width(Length::Fixed(220.0)), text_input("", value).on_input(move |s| Message::Settings(f(s))).padding(5).width(100)]
                .spacing(8)
                .align_y(iced::Center)
        };
        let body = column![
            text("Melting temperature (SantaLucia 1998 nearest neighbor)").size(15),
            field("Oligo concentration (nM)", &self.oligo, SettingsMsg::Oligo),
            field("Monovalent cations Na⁺/K⁺ (mM)", &self.na, SettingsMsg::Na),
            field("Divalent cations Mg²⁺ (mM)", &self.mg, SettingsMsg::Mg),
            field("dNTP (mM)", &self.dntp, SettingsMsg::Dntp),
            rule::horizontal(1),
            text("Inclusivity").size(15),
            field("3' end window (bases)", &self.window, SettingsMsg::Window),
            text("A mismatch within this many 3'-terminal bases counts as a 3' mismatch.").size(12).style(style::muted),
            rule::horizontal(1),
            text("Appearance").size(15),
            checkbox(self.dark).label("Dark theme").on_toggle(|b| Message::Settings(SettingsMsg::Dark(b))),
        ]
        .spacing(10);
        let buttons = row![
            button("Defaults").on_press(Message::Settings(SettingsMsg::Defaults)).style(button::text),
            Space::new().width(Fill),
            button("Cancel").on_press(Message::CloseModal).style(button::secondary),
            button("Save").on_press(Message::SaveSettings).style(button::primary),
        ]
        .spacing(10);
        dialog_frame("Settings", body.into(), buttons.into(), 480.0)
    }
}

// ---------------------------------------------------------------- common

/// Prominent "cannot be undone" notice for destructive confirmations.
fn irreversible<'a>(what: &'a str) -> Element<'a, Message> {
    container(
        column![
            text("This cannot be undone.").size(13).font(iced::Font { weight: iced::font::Weight::Bold, ..iced::Font::DEFAULT }),
            text(what).size(12),
        ]
        .spacing(2),
    )
    .padding([8, 10])
    .width(Fill)
    .style(|t: &iced::Theme| {
        let p = t.extended_palette();
        container::Style {
            background: Some(iced::Color { a: 0.12, ..p.danger.base.color }.into()),
            border: iced::Border { color: p.danger.base.color, width: 1.0, radius: 5.0.into() },
            text_color: Some(p.danger.strong.color),
            ..Default::default()
        }
    })
    .into()
}

pub fn dialog_frame<'a>(title: &'a str, body: Element<'a, Message>, buttons: Element<'a, Message>, width: f32) -> Element<'a, Message> {
    container(column![text(title).size(20), body, buttons].spacing(16))
        .padding(22)
        .width(Length::Fixed(width))
        .max_height(820.0)
        .style(style::card)
        .into()
}

pub fn view_modal(modal: &Modal) -> Element<'_, Message> {
    match modal {
        Modal::Sequence(d) => d.view(),
        Modal::Annotation(d) => d.view(),
        Modal::ConfirmDeleteAnnotation { id, name, detail } => {
            let body = column![
                text(format!("Delete the annotation “{name}”?")).size(14),
                text(detail.as_str()).size(12).style(style::muted),
                irreversible("The annotation is removed permanently."),
            ]
            .spacing(10);
            let buttons = row![
                Space::new().width(Fill),
                button("Cancel").on_press(Message::CloseModal).style(button::secondary),
                button("Delete").on_press(Message::DeleteAnnotation(*id)).style(button::danger),
            ]
            .spacing(10);
            dialog_frame("Delete annotation", body.into(), buttons.into(), 420.0)
        }
        Modal::Map(d) => d.view(),
        Modal::Export(d) => d.view(),
        Modal::Settings(d) => d.view(),
        Modal::Rename { text: t, folder, .. } => {
            let body = text_input("Name", t)
                .id(RENAME_INPUT)
                .on_input(|s| Message::Modal(ModalMsg::RenameInput(s)))
                .on_submit(Message::Modal(ModalMsg::RenameConfirm))
                .padding(8);
            let buttons = row![
                Space::new().width(Fill),
                button("Cancel").on_press(Message::CloseModal).style(button::secondary),
                button("Rename").on_press(Message::Modal(ModalMsg::RenameConfirm)).style(button::primary),
            ]
            .spacing(10);
            dialog_frame(if *folder { "Rename folder" } else { "Rename document" }, body.into(), buttons.into(), 420.0)
        }
        Modal::ConfirmDelete { id, name, detail } => {
            let body = column![
                text(format!("Delete “{name}” from the library?")).size(14),
                text(detail.as_str()).size(12).style(style::muted),
                irreversible("The document, its sequences and its annotations are removed permanently."),
            ]
            .spacing(10);
            let buttons = row![
                Space::new().width(Fill),
                button("Cancel").on_press(Message::CloseModal).style(button::secondary),
                button("Delete").on_press(Message::ConfirmDelete(*id)).style(button::danger),
            ]
            .spacing(10);
            dialog_frame("Delete document", body.into(), buttons.into(), 420.0)
        }
        Modal::ConfirmDeleteFolder { id, name, documents, subfolders } => {
            let contents = match (documents, subfolders) {
                (0, 0) => "The folder is empty.".to_string(),
                (d, s) => format!("Its {d} document(s) and {s} subfolder(s) are kept and move up one level."),
            };
            let body = column![
                text(format!("Delete the folder “{name}”?")).size(14),
                text(contents).size(12).style(style::muted),
                irreversible("The folder itself is removed permanently."),
            ]
            .spacing(10);
            let buttons = row![
                Space::new().width(Fill),
                button("Cancel").on_press(Message::CloseModal).style(button::secondary),
                button("Delete folder").on_press(Message::ConfirmDeleteFolder(*id)).style(button::danger),
            ]
            .spacing(10);
            dialog_frame("Delete folder", body.into(), buttons.into(), 440.0)
        }
        Modal::NameFields { delimiter, keys } => {
            let body = column![
                text("Splits every sequence name at a delimiter and stores the parts as metadata fields (for sorting and grouping).").size(13),
                row![
                    text("Delimiter").width(Length::Fixed(140.0)),
                    text_input("|", delimiter).on_input(|s| Message::Modal(ModalMsg::DelimiterInput(s))).padding(5).width(60)
                ]
                .align_y(iced::Center),
                row![
                    text("Field names").width(Length::Fixed(140.0)),
                    text_input("e.g. accession, genotype, country, year", keys).on_input(|s| Message::Modal(ModalMsg::KeysInput(s))).padding(5)
                ]
                .align_y(iced::Center),
                text("Comma separated; missing names become “Name field N”.").size(12).style(style::muted),
            ]
            .spacing(10);
            let buttons = row![
                Space::new().width(Fill),
                button("Cancel").on_press(Message::CloseModal).style(button::secondary),
                button("Apply").on_press(Message::Modal(ModalMsg::ApplyNameFields)).style(button::primary),
            ]
            .spacing(10);
            dialog_frame("Metadata from sequence names", body.into(), buttons.into(), 520.0)
        }
        Modal::Info { title, text: t } => {
            let body = text(t.as_str()).size(14);
            let buttons = row![Space::new().width(Fill), button("OK").on_press(Message::CloseModal).style(button::primary)];
            dialog_frame(title.as_str(), body.into(), buttons.into(), 480.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::state::ViewPrefs;
    use super::*;
    use crate::model::Document;

    fn open(rows: &[&str], reference: Option<usize>) -> OpenDoc {
        let doc = Document::new(rows.iter().enumerate().map(|(i, s)| Row::from_gapped(format!("s{i}"), String::new(), s.as_bytes().to_vec())).collect());
        let info = DocInfo {
            id: 1,
            name: "d".into(),
            kind: DocKind::Alignment,
            n_rows: doc.rows.len(),
            width: doc.width,
            reference,
            created: String::new(),
            description: String::new(),
            folder: None,
        };
        OpenDoc::new(info, doc, Vec::new(), &ViewPrefs::default())
    }

    #[test]
    fn annotation_range_fields() {
        // Reference positions; the reference has an insertion at columns 2-3.
        let od = open(&["AC--GTAC", "ACTTGTAC"], Some(0));
        let mut d = AnnotationDialog::new(None, "a".into(), AnnotationKind::Region, 1, 5, &od);
        assert_eq!((d.start.as_str(), d.end.as_str()), ("2", "3"));
        assert_eq!(d.range(), Ok((1, 5)));
        // A changed field is converted; the unchanged one keeps its column.
        d.end = "5".into();
        assert_eq!(d.range(), Ok((1, 7)));
        d.start = "4".into();
        d.end = "3".into();
        assert!(d.range().is_err());
        d.start = "1".into();
        d.end = "7".into();
        assert!(d.range().is_err());
        // Without a reference the fields are alignment columns.
        let od = open(&["ACGT", "ACGA"], None);
        let mut d = AnnotationDialog::new(Some(3), "b".into(), AnnotationKind::Probe, 0, 2, &od);
        assert_eq!((d.start.as_str(), d.end.as_str()), ("1", "2"));
        d.end = "4".into();
        assert_eq!(d.range(), Ok((0, 4)));
        d.name = " ".into();
        assert_eq!(d.problem().as_deref(), Some("Enter a name."));
    }
}
