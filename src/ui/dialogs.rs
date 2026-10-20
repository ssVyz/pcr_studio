//! Modal dialogs: new/edit sequence, map to reference, export, settings,
//! primersearch, text viewer, rename, confirmations.

use super::canvas::group_digits;
use super::job;
use super::state::OpenDoc;
use super::{AppSettings, Message, style};
use crate::fasta::{self, PastedRecord};
use crate::mapper::{MapParams, Sensitivity};
use crate::model::{AnnotationKind, DocInfo, DocKind, RefCoords, Row};
use crate::primer::Orientation;
use crate::seq::{self, GAP};
use crate::thermo::TmConditions;
use crate::tools::primersearch::{self as ps, Criterion, InputStats, Mode};
use crate::tools::{self, ToolConfig, ToolSource, ToolSpec, ToolStatus};
use iced::widget::{Space, button, checkbox, column, container, pick_list, row, rule, scrollable, text, text_editor, text_input};
use iced::{Element, Fill, Font, Length, Task};
use std::collections::HashMap;
use std::path::PathBuf;

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
    Primersearch(PrimersearchDialog),
    TextView(TextView),
    /// Exit requested while a tool runs; `detail` describes the run.
    ConfirmExit { detail: String },
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
            docs: library.iter().filter(|x| x.id != od.info.id && x.kind != DocKind::Text).map(|x| DocOption { id: x.id, name: x.name.clone() }).collect(),
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

/// Location setting of one external tool in the settings dialog.
pub struct ToolRow {
    pub spec: &'static ToolSpec,
    pub config: ToolConfig,
    /// Result of checking `config`; `None`: changed since the last check.
    pub status: Option<ToolStatus>,
}

pub struct SettingsDialog {
    oligo: String,
    na: String,
    mg: String,
    dntp: String,
    window: String,
    dark: bool,
    pub tools: Vec<ToolRow>,
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
    ToolSource(usize, ToolSource),
    ToolPath(usize, String),
    ToolBrowse(usize),
    ToolPicked(usize, Option<PathBuf>),
    ToolTest(usize),
    /// Result of a test of the given configuration.
    ToolTested(usize, ToolConfig, ToolStatus),
}

impl SettingsDialog {
    /// `status`: the current status of each tool (by id).
    pub fn new(s: &AppSettings, status: &HashMap<&'static str, ToolStatus>) -> SettingsDialog {
        let tm = s.primer.tm;
        SettingsDialog {
            oligo: format!("{}", tm.oligo_nm),
            na: format!("{}", tm.monovalent_mm),
            mg: format!("{}", tm.divalent_mm),
            dntp: format!("{}", tm.dntp_mm),
            window: s.primer.three_prime_window.to_string(),
            dark: s.dark,
            tools: tools::ALL
                .into_iter()
                .map(|spec| ToolRow { spec, config: s.tool_config(spec.id), status: status.get(spec.id).cloned() })
                .collect(),
        }
    }

    pub fn update(&mut self, m: SettingsMsg) -> Option<Task<Message>> {
        match m {
            SettingsMsg::Oligo(s) => self.oligo = s,
            SettingsMsg::Na(s) => self.na = s,
            SettingsMsg::Mg(s) => self.mg = s,
            SettingsMsg::Dntp(s) => self.dntp = s,
            SettingsMsg::Window(s) => self.window = s,
            SettingsMsg::Dark(b) => self.dark = b,
            SettingsMsg::Defaults => {
                // Tool locations are kept.
                let d = AppSettings { dark: self.dark, ..AppSettings::default() };
                let tools = std::mem::take(&mut self.tools);
                *self = SettingsDialog::new(&d, &HashMap::new());
                self.tools = tools;
            }
            SettingsMsg::ToolSource(i, src) => {
                if let Some(r) = self.tools.get_mut(i).filter(|r| r.config.source != src) {
                    r.config.source = src;
                    r.status = None;
                }
            }
            SettingsMsg::ToolPath(i, p) => {
                if let Some(r) = self.tools.get_mut(i) {
                    r.config.custom_path = p;
                    r.status = None;
                }
            }
            SettingsMsg::ToolBrowse(i) => {
                let name = self.tools.get(i)?.spec.name;
                return Some(Task::perform(
                    rfd::AsyncFileDialog::new().set_title(format!("Locate {name}")).add_filter("Programs", &["exe"]).add_filter("All files", &["*"]).pick_file(),
                    move |f| Message::Settings(SettingsMsg::ToolPicked(i, f.map(|f| f.path().to_path_buf()))),
                ));
            }
            SettingsMsg::ToolPicked(i, Some(path)) => {
                if let Some(r) = self.tools.get_mut(i) {
                    r.config = ToolConfig { source: ToolSource::Custom, custom_path: path.display().to_string() };
                    r.status = None;
                }
            }
            SettingsMsg::ToolPicked(_, None) => {}
            SettingsMsg::ToolTest(i) => {
                let r = self.tools.get_mut(i)?;
                r.status = Some(ToolStatus::Checking);
                let (spec, cfg) = (r.spec, r.config.clone());
                let tested = cfg.clone();
                return Some(job::background(
                    move || tools::check(spec, &tested),
                    move |st| {
                        let st = st.unwrap_or_else(|| ToolStatus::Unavailable("the check failed".into()));
                        Message::Settings(SettingsMsg::ToolTested(i, cfg, st))
                    },
                ));
            }
            SettingsMsg::ToolTested(i, cfg, st) => {
                if let Some(r) = self.tools.get_mut(i).filter(|r| r.config == cfg) {
                    r.status = Some(st);
                }
            }
        }
        None
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
        for r in &self.tools {
            s.tools.insert(r.spec.id.to_string(), r.config.clone());
        }
    }

    pub fn view(&self) -> Element<'_, Message> {
        let field = |label: &'static str, value: &str, f: fn(String) -> SettingsMsg| {
            row![text(label).width(Length::Fixed(220.0)), text_input("", value).on_input(move |s| Message::Settings(f(s))).padding(5).width(100)]
                .spacing(8)
                .align_y(iced::Center)
        };
        let mut body = column![
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
            rule::horizontal(1),
            text("External tools").size(15),
            text("Programs PCR Studio runs for some features. Without them those features are unavailable; everything else works.")
                .size(12)
                .style(style::muted),
        ]
        .spacing(10);
        for (i, r) in self.tools.iter().enumerate() {
            body = body.push(tool_row(i, r));
        }
        let buttons = row![
            button("Defaults").on_press(Message::Settings(SettingsMsg::Defaults)).style(button::text),
            Space::new().width(Fill),
            button("Cancel").on_press(Message::CloseModal).style(button::secondary),
            button("Save").on_press(Message::SaveSettings).style(button::primary),
        ]
        .spacing(10);
        let body = scrollable(body.padding(iced::Padding { right: 12.0, ..Default::default() })).height(Length::Shrink);
        dialog_frame("Settings", body.into(), buttons.into(), 580.0)
    }
}

fn tool_row(i: usize, r: &ToolRow) -> Element<'_, Message> {
    let set = Message::Settings;
    let bold = Font { weight: iced::font::Weight::Bold, ..Font::DEFAULT };
    let mut c = column![
        row![
            text(r.spec.name).size(14).font(bold).width(Length::Fixed(120.0)),
            pick_list(&ToolSource::ALL[..], Some(r.config.source), move |src| set(SettingsMsg::ToolSource(i, src))).text_size(13).width(Fill),
            button(text("Test").size(13)).padding([5, 12]).style(button::secondary).on_press(set(SettingsMsg::ToolTest(i))),
        ]
        .spacing(8)
        .align_y(iced::Center),
        text(r.spec.purpose).size(12).style(style::muted),
    ]
    .spacing(6);
    let exe = format!("{}{}", r.spec.exe, std::env::consts::EXE_SUFFIX);
    match r.config.source {
        ToolSource::Bundled => {
            c = c.push(text(format!("PCR Studio's own copy: {exe} in assets\\tools next to the program.")).size(12).style(style::muted));
        }
        ToolSource::Path => c = c.push(text(format!("The first {exe} found in the folders of the PATH variable.")).size(12).style(style::muted)),
        ToolSource::Custom => {
            c = c.push(
                row![
                    text_input("Path to the program", &r.config.custom_path).on_input(move |p| set(SettingsMsg::ToolPath(i, p))).padding(5).size(13),
                    button(text("Browse…").size(13)).padding([5, 10]).style(button::secondary).on_press(set(SettingsMsg::ToolBrowse(i))),
                ]
                .spacing(8)
                .align_y(iced::Center),
            );
        }
    }
    let status: Element<Message> = match &r.status {
        None => text("Not checked yet: press Test, or Save to check.").size(12).style(style::muted).into(),
        Some(ToolStatus::Checking) => text("Checking…").size(12).style(style::muted).into(),
        Some(ToolStatus::Ready { path, version }) => text(format!("✓ {} {version} · {}", r.spec.name, path.display()))
            .size(12)
            .wrapping(text::Wrapping::WordOrGlyph)
            .style(|t: &iced::Theme| text::Style { color: Some(t.extended_palette().success.base.color) })
            .into(),
        Some(ToolStatus::Unavailable(e)) => text(format!("✗ Not available: {e}"))
            .size(12)
            .wrapping(text::Wrapping::WordOrGlyph)
            .style(|t: &iced::Theme| text::Style { color: Some(t.extended_palette().danger.base.color) })
            .into(),
    };
    c.push(status).into()
}

// ---------------------------------------------------------------- primersearch

/// Numeric fields of the primersearch dialog.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PsField {
    Tm,
    ThreePrime,
    Target,
    MaxAmb,
    MaxSeeds,
    NOligos,
    Mismatches,
    Ambiguities,
    MaxCandidates,
    MaxWork,
    Threads,
    Oligo,
    Na,
    Mg,
    Dntp,
}

const PS_FIELDS: usize = 15;

/// Runs primersearch on the selected columns.
pub struct PrimersearchDialog {
    pub rows: ExportRows,
    pub name: String,
    pub mode: Mode,
    pub fixed: bool,
    pub orientation: Orientation,
    pub exclude_n: bool,
    pub only_twofold: bool,
    pub criterion: Criterion,
    pub spacer: bool,
    fields: [String; PS_FIELDS],
    /// Selected columns `[c0, c1)` and their description.
    pub c0: usize,
    pub c1: usize,
    pub range: String,
    /// The columns with gaps in the reference are left out (Collapse gaps).
    pub collapsed: bool,
    pub stats: InputStats,
    pub version: String,
    conditions: TmConditions,
    n_selected: usize,
    doc_name: String,
}

#[derive(Debug, Clone)]
pub enum PrimersearchMsg {
    Rows(ExportRows),
    Name(String),
    Mode(Mode),
    Fixed(bool),
    Orientation(Orientation),
    ExcludeN(bool),
    OnlyTwofold(bool),
    Criterion(Criterion),
    Spacer(bool),
    Field(PsField, String),
    Defaults,
}

impl PrimersearchDialog {
    /// `last`: the parameters of the previous run; `conditions`: the Tm settings.
    pub fn new(od: &OpenDoc, c0: usize, c1: usize, last: &ps::Params, conditions: TmConditions, version: &str) -> PrimersearchDialog {
        let range = match od.ref_range(c0, c1) {
            Some((a, b)) => format!(
                "ref {}–{} (columns {}–{})",
                group_digits(a as usize),
                group_digits(b as usize),
                group_digits(c0 + 1),
                group_digits(c1)
            ),
            None => format!("columns {}–{}", group_digits(c0 + 1), group_digits(c1)),
        };
        let short = match od.ref_range(c0, c1) {
            Some((a, b)) => format!("ref {a}-{b}"),
            None => format!("columns {}-{}", c0 + 1, c1),
        };
        let mut d = PrimersearchDialog {
            rows: if od.selected_rows.is_empty() { ExportRows::All } else { ExportRows::Selected },
            name: format!("{} [{short}] primersearch", od.info.name),
            mode: Mode::default(),
            fixed: false,
            orientation: Orientation::Forward,
            exclude_n: false,
            only_twofold: false,
            criterion: Criterion::default(),
            spacer: true,
            fields: Default::default(),
            c0,
            c1,
            range,
            collapsed: od.kept_columns().is_some(),
            stats: InputStats::default(),
            version: version.to_string(),
            conditions,
            n_selected: od.selected_rows.len(),
            doc_name: od.info.name.clone(),
        };
        d.set_params(last);
        d.refresh(od);
        d
    }

    fn set_params(&mut self, p: &ps::Params) {
        self.mode = p.mode;
        self.fixed = p.fixed;
        self.orientation = if p.reverse { Orientation::Reverse } else { Orientation::Forward };
        self.exclude_n = p.exclude_n;
        self.only_twofold = p.only_twofold;
        self.criterion = p.criterion;
        self.spacer = p.spacer;
        let c = &self.conditions;
        self.fields = [
            p.tm.to_string(),
            p.three_prime.to_string(),
            p.target.to_string(),
            p.max_amb.to_string(),
            p.max_seeds.to_string(),
            p.n_oligos.to_string(),
            p.mismatches.to_string(),
            p.ambiguities.to_string(),
            p.max_candidates.to_string(),
            p.max_work.to_string(),
            p.threads.to_string(),
            c.oligo_nm.to_string(),
            c.monovalent_mm.to_string(),
            c.divalent_mm.to_string(),
            c.dntp_mm.to_string(),
        ];
    }

    /// Recounts the input after the row choice changed.
    fn refresh(&mut self, od: &OpenDoc) {
        let rows = super::chosen_rows(od, self.rows);
        self.stats = ps::Input::new(&od.doc, &rows, self.c0, self.c1, od.kept_columns()).stats();
    }

    fn field(&self, f: PsField) -> &str {
        self.fields[f as usize].trim()
    }

    pub fn update(&mut self, m: PrimersearchMsg, od: Option<&OpenDoc>) {
        match m {
            PrimersearchMsg::Rows(r) => {
                self.rows = r;
                if let Some(od) = od {
                    self.refresh(od);
                }
            }
            PrimersearchMsg::Name(s) => self.name = s,
            PrimersearchMsg::Mode(x) => self.mode = x,
            PrimersearchMsg::Fixed(b) => self.fixed = b,
            PrimersearchMsg::Orientation(o) => self.orientation = o,
            PrimersearchMsg::ExcludeN(b) => self.exclude_n = b,
            PrimersearchMsg::OnlyTwofold(b) => self.only_twofold = b,
            PrimersearchMsg::Criterion(c) => self.criterion = c,
            PrimersearchMsg::Spacer(b) => self.spacer = b,
            PrimersearchMsg::Field(f, s) => self.fields[f as usize] = s,
            PrimersearchMsg::Defaults => self.set_params(&ps::Params::default()),
        }
    }

    /// The run parameters and conditions, or the first invalid field.
    pub fn params(&self) -> Result<(ps::Params, ps::Conditions), String> {
        let float = |f: PsField, what: &str, lo: f64, hi: f64| {
            self.field(f)
                .replace(',', ".")
                .parse::<f64>()
                .ok()
                .filter(|x| (lo..=hi).contains(x))
                .ok_or_else(|| format!("{what} must be a number from {lo} to {hi}."))
        };
        let int = |f: PsField, what: &str, min: u64| {
            self.field(f).replace([',', '_', ' '], "").parse::<u64>().ok().filter(|&x| x >= min).ok_or_else(|| {
                if min > 0 { format!("{what} must be a whole number of at least {min}.") } else { format!("{what} must be a whole number.") }
            })
        };
        let p = ps::Params {
            mode: self.mode,
            fixed: self.fixed,
            reverse: self.orientation == Orientation::Reverse,
            tm: float(PsField::Tm, "The minimum Tm", 0.0, 100.0)?,
            three_prime: int(PsField::ThreePrime, "Conserved 3' bases", 0)? as usize,
            target: float(PsField::Target, "The coverage target", 0.0, 100.0)?,
            max_amb: int(PsField::MaxAmb, "Ambiguity codes per primer", 0)? as usize,
            max_seeds: int(PsField::MaxSeeds, "Seeds per range", 0)? as usize,
            exclude_n: self.exclude_n,
            only_twofold: self.only_twofold,
            criterion: self.criterion,
            n_oligos: int(PsField::NOligos, "Oligos in the set", 1)? as usize,
            mismatches: int(PsField::Mismatches, "Mismatches", 0)? as usize,
            ambiguities: int(PsField::Ambiguities, "Ambiguity codes per oligo", 0)? as usize,
            max_candidates: int(PsField::MaxCandidates, "Max. candidates", 0)?,
            max_work: int(PsField::MaxWork, "Max. work", 0)?,
            spacer: self.spacer,
            threads: int(PsField::Threads, "Threads", 0)? as usize,
        };
        let c = ps::Conditions {
            oligo_um: float(PsField::Oligo, "The oligo concentration", 0.001, 1e9)? / 1000.0,
            na_mm: float(PsField::Na, "Na⁺", 0.0, 1e6)?,
            mg_mm: float(PsField::Mg, "Mg²⁺", 0.0, 1e6)?,
            dntp_mm: float(PsField::Dntp, "dNTP", 0.0, 1e6)?,
        };
        Ok((p, c))
    }

    /// Why the run cannot start yet.
    pub fn problem(&self) -> Option<String> {
        if self.stats.usable == 0 {
            return Some("No sequence is usable in the selected columns.".into());
        }
        if self.name.trim().is_empty() {
            return Some("Enter a name for the result.".into());
        }
        self.params().err()
    }

    pub fn view(&self) -> Element<'_, Message> {
        let m = Message::Primersearch;
        let muted = |t: String| text(t).size(12).style(style::muted).wrapping(text::Wrapping::WordOrGlyph);
        let field = |label: &'static str, f: PsField| {
            row![
                text(label).size(13).width(Fill),
                text_input("", &self.fields[f as usize]).on_input(move |s| m(PrimersearchMsg::Field(f, s))).padding(4).size(13).width(120)
            ]
            .spacing(8)
            .align_y(iced::Center)
        };
        let labeled = |label: &'static str, w: Element<'static, Message>| {
            row![text(label).size(13).width(Length::Fixed(90.0)), w].spacing(8).align_y(iced::Center)
        };
        let s = &self.stats;
        let d = group_digits;
        let mut options = vec![ExportRows::All, ExportRows::Shown];
        if self.n_selected > 0 {
            options.push(ExportRows::Selected);
        }
        let empty = s.usable == 0;
        let usable = text(format!("{} of {} sequences usable", d(s.usable), d(s.written)))
            .size(13)
            .font(Font { weight: iced::font::Weight::Semibold, ..Font::DEFAULT })
            .style(move |t: &iced::Theme| text::Style { color: empty.then(|| t.extended_palette().danger.base.color) });
        let mut dropped = format!(
            "primersearch drops sequences with a gap, missing coverage or an N/ambiguity code in the selection: \
             {} with gaps or no coverage, {} with N/ambiguity codes.",
            d(s.gaps),
            d(s.ambiguous)
        );
        if s.chosen > s.written {
            dropped.push_str(&format!(" {} without any base in the selection are left out.", d(s.chosen - s.written)));
        }
        let mut input = column![
            text("Input").size(15),
            text(format!("Selection: {}, {} columns", self.range, d(s.columns))).size(13).wrapping(text::Wrapping::WordOrGlyph),
        ]
        .spacing(8);
        if self.collapsed {
            input = input.push(muted("Columns with gaps in the reference are left out (Collapse gaps).".into()));
        }
        input = input.push(pick_list(options, Some(self.rows), move |r| m(PrimersearchMsg::Rows(r))).text_size(13).width(Fill));
        input = input.push(usable).push(muted(dropped));
        let left = column![
            input,
            rule::horizontal(1),
            text("Reaction conditions").size(15),
            field("Oligo concentration (nM)", PsField::Oligo),
            field("Na⁺/K⁺ (mM)", PsField::Na),
            field("Mg²⁺ (mM)", PsField::Mg),
            field("dNTP (mM)", PsField::Dntp),
            muted("From Settings; changes here apply to this run only.".into()),
            rule::horizontal(1),
            text("Result").size(15),
            text_input("Name of the result", &self.name).on_input(move |s| m(PrimersearchMsg::Name(s))).padding(6),
            muted(format!(
                "Saved as a text document in the library (in the folder of “{}”) when the run finishes. \
                 PCR Studio stays usable while it runs.",
                self.doc_name
            )),
        ]
        .spacing(8)
        .width(Fill);

        let mut right = column![
            text("Search").size(15),
            labeled("Mode", pick_list(&Mode::ALL[..], Some(self.mode), move |x| m(PrimersearchMsg::Mode(x))).text_size(13).width(Fill).into()),
            checkbox(self.fixed).label("Fixed slice: the selection is the primer site").on_toggle(move |b| m(PrimersearchMsg::Fixed(b))).text_size(13),
            muted(if self.fixed {
                "Generates the primer variants that cover all sequences at the selected columns. The Tm threshold is not enforced.".into()
            } else {
                "Searches the selected columns for primers that reach the minimum Tm.".into()
            }),
            labeled(
                "Orientation",
                pick_list([Orientation::Forward, Orientation::Reverse], Some(self.orientation), move |o| m(PrimersearchMsg::Orientation(o)))
                    .text_size(13)
                    .width(Fill)
                    .into()
            ),
            field("Minimum Tm (°C)", PsField::Tm),
            field("Conserved 3' bases", PsField::ThreePrime),
        ]
        .spacing(8)
        .width(Fill);
        match self.mode {
            Mode::NoAmbiguities => {}
            Mode::Incremental => {
                right = right
                    .push(field("Coverage target (%)", PsField::Target))
                    .push(field("Ambiguity codes per primer", PsField::MaxAmb))
                    .push(field("Seeds per range (0 = no cap)", PsField::MaxSeeds));
            }
            Mode::OptimizeByMismatch => {
                right = right
                    .push(
                        row![
                            text("Count sequences with").size(13).width(Fill),
                            pick_list(&Criterion::ALL[..], Some(self.criterion), move |c| m(PrimersearchMsg::Criterion(c))).text_size(13).width(200)
                        ]
                        .spacing(8)
                        .align_y(iced::Center),
                    )
                    .push(field("Mismatches", PsField::Mismatches))
                    .push(field("Oligos in the set", PsField::NOligos))
                    .push(field("Ambiguity codes per oligo", PsField::Ambiguities))
                    .push(field("Max. candidates (0 = no limit)", PsField::MaxCandidates))
                    .push(field("Max. work (0 = no limit)", PsField::MaxWork))
                    .push(muted("The search is exhaustive; a run that exceeds a limit stops with an explanation.".into()));
            }
        }
        if self.mode.uses_iupac() {
            right = right
                .push(checkbox(self.exclude_n).label("No N in primers").on_toggle(move |b| m(PrimersearchMsg::ExcludeN(b))).text_size(13))
                .push(
                    checkbox(self.only_twofold)
                        .label("Only two-fold codes (R, Y, S, W, K, M)")
                        .on_toggle(move |b| m(PrimersearchMsg::OnlyTwofold(b)))
                        .text_size(13),
                );
        }
        right = right
            .push(rule::horizontal(1))
            .push(text("Output").size(15))
            .push(checkbox(self.spacer).label("Primer sequences in groups of three bases").on_toggle(move |b| m(PrimersearchMsg::Spacer(b))).text_size(13))
            .push(field("Threads (0 = all cores)", PsField::Threads));

        let body = scrollable(row![left, right].spacing(32).padding(iced::Padding { right: 12.0, ..Default::default() })).height(Length::Shrink);
        let problem = self.problem();
        let buttons = row![
            text(problem.clone().unwrap_or_default()).size(12).style(style::muted).width(Fill),
            button("Defaults").on_press(m(PrimersearchMsg::Defaults)).style(button::text),
            button("Cancel").on_press(Message::CloseModal).style(button::secondary),
            button("Run").on_press_maybe(problem.is_none().then_some(Message::RunPrimersearch)).style(button::primary),
        ]
        .spacing(10)
        .align_y(iced::Center);
        dialog_frame_owned(format!("Find primers · primersearch {}", self.version), body.into(), buttons.into(), 920.0)
    }
}

// ---------------------------------------------------------------- text viewer

/// Shows a text document (read-only; text can be selected and copied).
pub struct TextView {
    pub name: String,
    pub info: String,
    pub content: text_editor::Content,
}

#[derive(Debug, Clone)]
pub enum TextViewMsg {
    Action(text_editor::Action),
    Copy,
    SaveAs,
}

impl TextView {
    pub fn new(name: String, info: String, content: &str) -> TextView {
        TextView { name, info, content: text_editor::Content::with_text(content) }
    }

    pub fn update(&mut self, m: TextViewMsg) -> Option<Task<Message>> {
        match m {
            // Selecting and scrolling only; the text cannot be changed.
            TextViewMsg::Action(a) => {
                if !a.is_edit() {
                    self.content.perform(a);
                }
                None
            }
            TextViewMsg::Copy => Some(iced::clipboard::write(self.content.text())),
            TextViewMsg::SaveAs => Some(Task::perform(
                rfd::AsyncFileDialog::new()
                    .set_title("Save as text file")
                    .set_file_name(format!("{}.txt", sanitize(&self.name)))
                    .add_filter("Text", &["txt"])
                    .save_file(),
                |f| Message::TextSavePicked(f.map(|f| f.path().to_path_buf())),
            )),
        }
    }

    pub fn view(&self) -> Element<'_, Message> {
        let t = Message::TextView;
        let body = column![
            text(self.info.as_str()).size(12).style(style::muted).wrapping(text::Wrapping::WordOrGlyph),
            text_editor(&self.content)
                .on_action(move |a| t(TextViewMsg::Action(a)))
                .font(Font::MONOSPACE)
                .size(12)
                .padding(8)
                .height(Length::Fixed(600.0))
                .wrapping(text::Wrapping::WordOrGlyph),
        ]
        .spacing(10);
        let buttons = row![
            button("Copy all").on_press(t(TextViewMsg::Copy)).style(button::secondary),
            button("Save as .txt…").on_press(t(TextViewMsg::SaveAs)).style(button::secondary),
            Space::new().width(Fill),
            button("Close").on_press(Message::CloseModal).style(button::primary),
        ]
        .spacing(10);
        dialog_frame(self.name.as_str(), body.into(), buttons.into(), 1040.0)
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
    dialog_frame_owned(title.to_string(), body, buttons, width)
}

pub fn dialog_frame_owned<'a>(title: String, body: Element<'a, Message>, buttons: Element<'a, Message>, width: f32) -> Element<'a, Message> {
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
            let body = text(t.as_str()).size(14).wrapping(text::Wrapping::WordOrGlyph);
            let buttons = row![Space::new().width(Fill), button("OK").on_press(Message::CloseModal).style(button::primary)];
            dialog_frame(title.as_str(), body.into(), buttons.into(), 480.0)
        }
        Modal::Primersearch(d) => d.view(),
        Modal::TextView(v) => v.view(),
        Modal::ConfirmExit { detail } => {
            let body = column![
                text(detail.as_str()).size(14),
                text("Exiting stops it; its output is not saved.").size(13).style(style::muted),
            ]
            .spacing(10);
            let buttons = row![
                Space::new().width(Fill),
                button("Keep running").on_press(Message::CloseModal).style(button::secondary),
                button("Stop and exit").on_press(Message::ConfirmExit).style(button::danger),
            ]
            .spacing(10);
            dialog_frame("Exit PCR Studio?", body.into(), buttons.into(), 460.0)
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
