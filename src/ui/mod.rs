//! The iced application: state, messages and update logic.

pub mod canvas;
mod icons;
mod menu;
mod dialogs;
mod job;
pub mod state;
pub mod style;
mod view;

use crate::db::{self, Db};
use crate::display::SortKey;
use crate::fasta;
use crate::mapper::{self, MapParams};
use crate::model::{Annotation, AnnotationKind, ConsensusThreshold, DocInfo, DocKind, Document, Folder, Row};
use crate::primer::{self, OligoSource, Orientation, PrimerOptions};
use crate::search::{self, Scope};
use canvas::ViewMsg;
use dialogs::{ExportDialog, ExportRows, MapDialog, Modal, RefChoice, SettingsDialog};
use iced::keyboard::{self, Key, key::Named};
use iced::{Subscription, Task, Theme};
use job::{Job, Payload};
use serde::{Deserialize, Serialize};
use state::{BASE_COL_W, Highlight, OpenDoc, Selection, ViewPrefs};
use std::path::PathBuf;
use std::sync::Arc;

pub fn run() -> iced::Result {
    iced::application(App::boot, App::update, App::view)
        .title(App::title)
        .theme(App::theme)
        .subscription(App::subscription)
        .window_size((1500.0, 920.0))
        .antialiasing(false)
        .run()
}

/// Settings persisted in the database.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct AppSettings {
    pub view: ViewPrefs,
    pub primer: PrimerOptions,
    pub map: MapParams,
    pub dark: bool,
    pub row_height: f32,
}

impl Default for AppSettings {
    fn default() -> Self {
        AppSettings { view: ViewPrefs::default(), primer: PrimerOptions::default(), map: MapParams::default(), dark: false, row_height: 16.0 }
    }
}

const SETTINGS_KEY: &str = "app_settings";

/// Ids of the viewer's text fields (focused from the menus / shortcuts).
pub const JUMP_INPUT: &str = "jump-input";
pub const SEARCH_INPUT: &str = "search-input";

/// "Group by" choice for the pick list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroupChoice {
    None,
    Key(String),
}

impl std::fmt::Display for GroupChoice {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GroupChoice::None => f.write_str("No grouping"),
            GroupChoice::Key(k) => f.write_str(k),
        }
    }
}

/// A library entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LibKey {
    Doc(i64),
    Folder(i64),
}

/// Where a dragged library entry would land.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropTarget {
    Folder(i64),
    Root,
}

/// A press on a library entry; it becomes a drag once the cursor moves a few pixels.
#[derive(Debug, Clone)]
pub struct LibDrag {
    pub item: LibKey,
    origin: Option<iced::Point>,
    pub pos: Option<iced::Point>,
    pub active: bool,
    /// Collapsed folder being hovered during the drag, and since when (auto-expand).
    hover_since: Option<(i64, std::time::Instant)>,
    /// -1 / +1 while the cursor is in the top / bottom auto-scroll zone of the list.
    zone: i8,
}

const DRAG_THRESHOLD: f32 = 6.0;
const AUTO_EXPAND: std::time::Duration = std::time::Duration::from_millis(650);
pub const LIB_SCROLL: &str = "library-scroll";

/// Destination folder in the "Move to" pick list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FolderChoice {
    pub id: Option<i64>,
    pub label: String,
}

impl std::fmt::Display for FolderChoice {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.label)
    }
}

pub struct LoadedDoc {
    open: OpenDoc,
}

#[derive(Debug, Clone)]
pub enum Message {
    Tick,
    Key(keyboard::Event),
    // Library
    OpenDocument(i64),
    DocLoaded(Payload<Result<LoadedDoc, String>>),
    ImportFasta,
    ImportPicked(Option<Vec<PathBuf>>),
    ImportDone(Payload<Result<Vec<i64>, String>>),
    AskRename(i64),
    AskDelete(i64),
    ConfirmDelete(i64),
    ToggleFolder(i64),
    /// Press on a library entry: select it; may start a drag.
    LibPress(LibKey),
    /// Double-click: open a document / expand or collapse a folder.
    LibActivate(LibKey),
    LibHover(LibKey),
    LibUnhover(LibKey),
    DragMove(iced::Point),
    DragRelease,
    DragTick,
    DragZone(i8),
    NewFolder,
    AskRenameFolder(i64),
    AskDeleteFolder(i64),
    ConfirmDeleteFolder(i64),
    MoveTo(FolderChoice),
    ExtractSlice,
    SliceDone(Payload<Result<(i64, String), String>>),
    // Viewer
    View(ViewMsg),
    ZoomIn,
    ZoomOut,
    ZoomFit,
    ZoomLetters,
    RowHeight(f32),
    JumpInput(String),
    JumpSubmit,
    FilterInput(String),
    SearchInput(String),
    SearchSubmit,
    SearchStep(i32),
    SearchScope(Scope),
    SearchMismatches(u8),
    // Display preferences
    SetHighlight(Highlight),
    ToggleDots(bool),
    ToggleGaps(bool),
    ToggleConsensus(bool),
    ToggleIdentity(bool),
    ToggleAnnotations(bool),
    ToggleOverview(bool),
    ToggleColors(bool),
    SetGraph(state::GraphMode),
    SetThreshold(ConsensusThreshold),
    // Rows
    SortBy(SortKey),
    SortDescending(bool),
    GroupBy(GroupChoice),
    Collapse(bool),
    SetReferenceFromSelection,
    ClearReference,
    ClearRowSelection,
    // Selection / primer popup
    ClearSelection,
    ShowPopup,
    ClosePopup,
    PrimerSource(OligoSource),
    PrimerOrientation(Orientation),
    CopyOligo,
    CopyReport,
    AnnotationName(String),
    AnnotationKind(AnnotationKind),
    AddAnnotation,
    DeleteAnnotation(i64),
    GotoAnnotation(usize),
    // Metadata
    ImportMetadata,
    MetadataPicked(Option<PathBuf>),
    OpenNameFields,
    // Export
    OpenExport,
    Export(dialogs::ExportMsg),
    ExportPicked(Option<PathBuf>),
    ExportDone(Payload<Result<String, String>>),
    // Mapping
    OpenMapDialog,
    Map(dialogs::MapMsg),
    RunMapping,
    MapDone(Payload<Result<(i64, String), String>>),
    // Settings / modals
    OpenSettings,
    Settings(dialogs::SettingsMsg),
    SaveSettings,
    Modal(dialogs::ModalMsg),
    CloseModal,
    CancelJob,
    DemoData,
    // Menu bar
    ToggleMenu(menu::MenuId),
    HoverMenu(menu::MenuId),
    CloseMenu,
    /// A menu entry was chosen: close the menu, then run the action.
    MenuAction(Box<Message>),
    CloseDocument,
    FocusJump,
    FocusSearch,
    Exit,
    Noop,
}

pub struct App {
    db: Db,
    db_path: PathBuf,
    library: Vec<DocInfo>,
    library_selected: Option<i64>,
    folders: Vec<Folder>,
    /// Selected folder (mutually exclusive with `library_selected`).
    folder_selected: Option<i64>,
    /// Open drop-down menu.
    menu: Option<menu::MenuId>,
    /// Pending or active drag of a library entry.
    drag: Option<LibDrag>,
    /// Library entry under the cursor.
    lib_hover: Option<LibKey>,
    open: Option<OpenDoc>,
    settings: AppSettings,
    job: Option<Job>,
    modal: Option<Modal>,
    status: String,
    error: Option<String>,
}

impl App {
    fn boot() -> (App, Task<Message>) {
        let db_path = db::default_path();
        let (db, error) = match Db::open(&db_path) {
            Ok(db) => (db, None),
            Err(e) => {
                let fallback = std::env::temp_dir().join(db::DB_FILE_NAME);
                let db = Db::open(&fallback).expect("cannot open any database");
                (db, Some(format!("{e}. Using {} instead.", fallback.display())))
            }
        };
        let settings: AppSettings = db.setting(SETTINGS_KEY).and_then(|s| serde_json::from_str(&s).ok()).unwrap_or_default();
        let library = db.list_documents().unwrap_or_default();
        let folders = db.list_folders().unwrap_or_default();
        let status = format!("Library: {} · {} document(s)", db_path.display(), library.len());
        // FASTA files given on the command line ("Open with…") are imported right away.
        let paths: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).filter(|p| p.is_file()).collect();
        let task = if paths.is_empty() { Task::none() } else { Task::done(Message::ImportPicked(Some(paths))) };
        let app = App {
            db,
            db_path,
            library,
            library_selected: None,
            folders,
            folder_selected: None,
            menu: None,
            drag: None,
            lib_hover: None,
            open: None,
            settings,
            job: None,
            modal: None,
            status,
            error,
        };
        (app, task)
    }

    fn title(&self) -> String {
        match &self.open {
            Some(od) => format!("{} — PCR Studio", od.info.name),
            None => "PCR Studio".into(),
        }
    }

    fn theme(&self) -> Theme {
        style::app_theme(self.settings.dark)
    }

    fn subscription(&self) -> Subscription<Message> {
        let mut subs = vec![keyboard::listen().map(Message::Key)];
        if self.job.is_some() {
            subs.push(iced::time::every(std::time::Duration::from_millis(120)).map(|_| Message::Tick));
        }
        if let Some(d) = &self.drag {
            // Follow the cursor anywhere in the window until the button is released.
            subs.push(iced::event::listen_with(|event, _status, _window| match event {
                iced::Event::Mouse(iced::mouse::Event::CursorMoved { position }) => Some(Message::DragMove(position)),
                iced::Event::Mouse(iced::mouse::Event::ButtonReleased(iced::mouse::Button::Left)) => Some(Message::DragRelease),
                _ => None,
            }));
            if d.active {
                subs.push(iced::time::every(std::time::Duration::from_millis(60)).map(|_| Message::DragTick));
            }
        }
        Subscription::batch(subs)
    }

    /// Target under the cursor for the active drag (`None`: outside the library).
    fn drop_target(&self) -> Option<DropTarget> {
        let d = self.drag.as_ref().filter(|d| d.active)?;
        match self.lib_hover {
            Some(LibKey::Folder(f)) => Some(DropTarget::Folder(f)),
            Some(LibKey::Doc(id)) => {
                let folder = self.library.iter().find(|x| x.id == id).and_then(|x| x.folder);
                Some(folder.map_or(DropTarget::Root, DropTarget::Folder))
            }
            None => {
                let p = d.pos?;
                (p.x < view::LIBRARY_W && p.y > menu::BAR_H).then_some(DropTarget::Root)
            }
        }
    }

    /// Whether dropping `item` on `target` would move it.
    fn drop_allowed(&self, item: LibKey, target: DropTarget) -> bool {
        let dest = match target {
            DropTarget::Folder(f) => Some(f),
            DropTarget::Root => None,
        };
        match item {
            LibKey::Doc(id) => self.library.iter().find(|x| x.id == id).is_some_and(|x| x.folder != dest),
            LibKey::Folder(id) => {
                let parent = self.folders.iter().find(|f| f.id == id).and_then(|f| f.parent);
                parent != dest && dest.is_none_or(|f| !self.folder_subtree(id).contains(&f))
            }
        }
    }

    fn lib_name(&self, key: LibKey) -> String {
        match key {
            LibKey::Doc(id) => self.library.iter().find(|d| d.id == id).map(|d| d.name.clone()),
            LibKey::Folder(id) => self.folders.iter().find(|f| f.id == id).map(|f| f.name.clone()),
        }
        .unwrap_or_default()
    }

    fn finish_drag(&mut self) {
        let target = self.drop_target();
        let Some(d) = self.drag.take() else { return };
        if !d.active {
            return;
        }
        let Some(target) = target.filter(|&t| self.drop_allowed(d.item, t)) else {
            return;
        };
        let dest = match target {
            DropTarget::Folder(f) => Some(f),
            DropTarget::Root => None,
        };
        let result = match d.item {
            LibKey::Doc(id) => self.db.move_document(id, dest),
            LibKey::Folder(id) => self.db.move_folder(id, dest),
        };
        match result {
            Ok(()) => {
                let name = self.lib_name(d.item);
                let place = dest.map(|f| self.lib_name(LibKey::Folder(f))).unwrap_or_else(|| "the top level".into());
                self.reload_library();
                self.reveal_folder(dest);
                self.status = format!("Moved “{name}” to {place}");
            }
            Err(e) => self.error = Some(e),
        }
    }

    fn palette(&self) -> style::Palette {
        if self.settings.dark { style::Palette::dark() } else { style::Palette::light() }
    }

    fn save_settings(&mut self) {
        if let Ok(s) = serde_json::to_string(&self.settings)
            && let Err(e) = self.db.set_setting(SETTINGS_KEY, &s) {
                self.error = Some(e);
            }
    }

    fn reload_library(&mut self) {
        match self.db.list_documents() {
            Ok(l) => self.library = l,
            Err(e) => self.error = Some(e),
        }
        match self.db.list_folders() {
            Ok(f) => self.folders = f,
            Err(e) => self.error = Some(e),
        }
        if self.folder_selected.is_some_and(|id| !self.folders.iter().any(|f| f.id == id)) {
            self.folder_selected = None;
        }
    }

    /// Folder that receives new documents: the selected folder, or the folder
    /// of the selected document.
    fn target_folder(&self) -> Option<i64> {
        self.folder_selected
            .or_else(|| self.library_selected.and_then(|id| self.library.iter().find(|d| d.id == id)).and_then(|d| d.folder))
    }

    /// Makes all ancestors of `folder` expanded so its content is visible.
    fn reveal_folder(&mut self, folder: Option<i64>) {
        let mut cur = folder;
        while let Some(id) = cur {
            let Some(f) = self.folders.iter_mut().find(|f| f.id == id) else { break };
            if f.collapsed {
                f.collapsed = false;
                let _ = self.db.set_folder_collapsed(id, false);
            }
            cur = f.parent;
        }
    }

    /// Folder ids of `id` and everything below it.
    fn folder_subtree(&self, id: i64) -> Vec<i64> {
        let mut out = vec![id];
        let mut i = 0;
        while i < out.len() {
            let cur = out[i];
            out.extend(self.folders.iter().filter(|f| f.parent == Some(cur)).map(|f| f.id));
            i += 1;
        }
        out
    }

    fn busy(&self) -> bool {
        self.job.is_some()
    }

    fn start_job<T: Send + 'static>(
        &mut self,
        title: &str,
        cancellable: bool,
        work: impl FnOnce(&job::Progress) -> T + Send + 'static,
        done: impl FnOnce(Option<T>) -> Message + Send + 'static,
    ) -> Task<Message> {
        let (j, task) = job::spawn(title, cancellable, work, done);
        self.job = Some(j);
        task
    }

    fn open_document(&mut self, id: i64) -> Task<Message> {
        if self.busy() {
            return Task::none();
        }
        self.library_selected = Some(id);
        let path = self.db_path.clone();
        let prefs = self.settings.view.clone();
        let row_h = self.settings.row_height;
        self.start_job(
            "Opening document",
            false,
            move |p| {
                let db = Db::open(&path)?;
                let (info, doc, anns) = db.load_document(id, &|f| p.set(f * 0.8, "Loading sequences"))?;
                p.set(0.85, "Computing consensus and statistics");
                let mut open = OpenDoc::new(info, doc, anns, &prefs);
                open.set_row_h(row_h);
                Ok(LoadedDoc { open })
            },
            |r| Message::DocLoaded(Payload::new(r.unwrap_or_else(|| Err("Opening the document failed".into())))),
        )
    }

    fn with_doc(&mut self, f: impl FnOnce(&mut OpenDoc, &AppSettings)) {
        if let Some(od) = self.open.as_mut() {
            f(od, &self.settings);
        }
    }

    fn after_selection_change(&mut self) {
        let opts = self.settings.primer;
        if let Some(od) = self.open.as_mut() {
            od.evaluate_selection(&opts);
        }
    }

    fn update(&mut self, message: Message) -> Task<Message> {
        match message {
            Message::Tick | Message::Noop => {}
            Message::Key(e) => return self.on_key(e),
            Message::LibPress(key) => {
                match key {
                    LibKey::Doc(id) => {
                        self.library_selected = Some(id);
                        self.folder_selected = None;
                    }
                    LibKey::Folder(id) => {
                        self.folder_selected = Some(id);
                        self.library_selected = None;
                    }
                }
                self.drag = Some(LibDrag { item: key, origin: None, pos: None, active: false, hover_since: None, zone: 0 });
            }
            Message::LibActivate(key) => {
                self.drag = None;
                match key {
                    LibKey::Doc(id) => return self.open_document(id),
                    LibKey::Folder(id) => return self.update(Message::ToggleFolder(id)),
                }
            }
            Message::LibHover(key) => {
                self.lib_hover = Some(key);
                if let (Some(d), LibKey::Folder(f)) = (self.drag.as_mut(), key) {
                    d.hover_since = Some((f, std::time::Instant::now()));
                }
            }
            Message::LibUnhover(key) => {
                if self.lib_hover == Some(key) {
                    self.lib_hover = None;
                }
                if let (Some(d), LibKey::Folder(f)) = (self.drag.as_mut(), key)
                    && d.hover_since.is_some_and(|(h, _)| h == f)
                {
                    d.hover_since = None;
                }
            }
            Message::DragMove(p) => {
                if let Some(d) = self.drag.as_mut() {
                    let origin = *d.origin.get_or_insert(p);
                    d.pos = Some(p);
                    if !d.active && ((p.x - origin.x).powi(2) + (p.y - origin.y).powi(2)).sqrt() > DRAG_THRESHOLD {
                        d.active = true;
                    }
                }
            }
            Message::DragRelease => self.finish_drag(),
            Message::DragZone(z) => {
                if let Some(d) = self.drag.as_mut() {
                    d.zone = z;
                }
            }
            Message::DragTick => {
                let Some(d) = self.drag.as_mut().filter(|d| d.active) else { return Task::none() };
                // Hovering a collapsed folder expands it, so items can be dropped deeper.
                if let Some((f, since)) = d.hover_since
                    && since.elapsed() >= AUTO_EXPAND
                {
                    d.hover_since = None;
                    if let Some(folder) = self.folders.iter_mut().find(|x| x.id == f && x.collapsed) {
                        folder.collapsed = false;
                        let _ = self.db.set_folder_collapsed(f, false);
                    }
                }
                let zone = self.drag.as_ref().map(|d| d.zone).unwrap_or(0);
                if zone != 0 {
                    return iced::widget::operation::scroll_by(LIB_SCROLL, iced::widget::scrollable::AbsoluteOffset { x: 0.0, y: zone as f32 * 18.0 });
                }
            }
            Message::ToggleFolder(id) => {
                if let Some(f) = self.folders.iter_mut().find(|f| f.id == id) {
                    f.collapsed = !f.collapsed;
                    if let Err(e) = self.db.set_folder_collapsed(id, f.collapsed) {
                        self.error = Some(e);
                    }
                }
            }
            Message::NewFolder => {
                let parent = self.target_folder();
                match self.db.create_folder("New folder", parent) {
                    Ok(id) => {
                        self.reload_library();
                        self.reveal_folder(parent);
                        self.folder_selected = Some(id);
                        self.library_selected = None;
                        self.modal = Some(Modal::Rename { id, text: "New folder".into(), folder: true });
                        return dialogs::focus_rename();
                    }
                    Err(e) => self.error = Some(e),
                }
            }
            Message::AskRenameFolder(id) => {
                if let Some(f) = self.folders.iter().find(|f| f.id == id) {
                    self.modal = Some(Modal::Rename { id, text: f.name.clone(), folder: true });
                    return dialogs::focus_rename();
                }
            }
            Message::AskDeleteFolder(id) => {
                if let Some(f) = self.folders.iter().find(|f| f.id == id) {
                    self.modal = Some(Modal::ConfirmDeleteFolder { id, name: f.name.clone() });
                }
            }
            Message::ConfirmDeleteFolder(id) => {
                self.modal = None;
                if let Err(e) = self.db.delete_folder(id) {
                    self.error = Some(e);
                }
                self.folder_selected = None;
                self.reload_library();
            }
            Message::MoveTo(choice) => {
                let result = match (self.library_selected, self.folder_selected) {
                    (Some(doc), _) => self.db.move_document(doc, choice.id),
                    (None, Some(folder)) => self.db.move_folder(folder, choice.id),
                    _ => Ok(()),
                };
                match result {
                    Ok(()) => {
                        self.reload_library();
                        self.reveal_folder(choice.id);
                        self.status = format!("Moved to {}", choice.label.trim_start_matches(['·', ' ']));
                    }
                    Err(e) => self.error = Some(e),
                }
            }
            Message::ExtractSlice => return self.extract_slice(),
            Message::SliceDone(p) => {
                self.job = None;
                match p.take() {
                    Some(Ok((id, msg))) => {
                        self.reload_library();
                        let folder = self.library.iter().find(|d| d.id == id).and_then(|d| d.folder);
                        self.reveal_folder(folder);
                        self.library_selected = Some(id);
                        self.folder_selected = None;
                        self.status = msg;
                    }
                    Some(Err(e)) => self.error = Some(e),
                    None => {}
                }
            }
            Message::OpenDocument(id) => return self.open_document(id),
            Message::DocLoaded(p) => {
                self.job = None;
                match p.take() {
                    Some(Ok(mut loaded)) => {
                        if let Some(old) = &self.open {
                            loaded.open.vp.canvas = old.vp.canvas;
                            loaded.open.vp.seq_w = old.vp.seq_w;
                            loaded.open.vp.rows_h = old.vp.rows_h;
                        }
                        let od = loaded.open;
                        self.status = format!(
                            "Opened “{}”: {} sequences, {} columns",
                            od.info.name,
                            canvas::group_digits(od.doc.rows.len()),
                            canvas::group_digits(od.doc.width)
                        );
                        self.library_selected = Some(od.info.id);
                        self.open = Some(od);
                    }
                    Some(Err(e)) => self.error = Some(e),
                    None => {}
                }
            }
            Message::ImportFasta => {
                if self.busy() {
                    return Task::none();
                }
                return Task::perform(
                    rfd::AsyncFileDialog::new()
                        .set_title("Import FASTA")
                        .add_filter("FASTA", &["fasta", "fa", "fas", "fna", "ffn", "faa", "aln", "mfa", "txt"])
                        .add_filter("All files", &["*"])
                        .pick_files(),
                    |files| Message::ImportPicked(files.map(|v| v.into_iter().map(|f| f.path().to_path_buf()).collect())),
                );
            }
            Message::ImportPicked(Some(paths)) if !paths.is_empty() => {
                let db_path = self.db_path.clone();
                let folder = self.target_folder();
                return self.start_job(
                    "Importing FASTA",
                    false,
                    move |p| {
                        let mut db = Db::open(&db_path)?;
                        let mut ids = Vec::new();
                        let n = paths.len() as f32;
                        for (i, path) in paths.iter().enumerate() {
                            let base = i as f32 / n;
                            let label = format!("Reading {}", path.file_name().map(|s| s.to_string_lossy()).unwrap_or_default());
                            let imp = fasta::read_fasta(path, &|f| p.set(base + f * 0.7 / n, &label))?;
                            p.set(base + 0.7 / n, "Saving to library");
                            let desc = format!("Imported from {}", path.display());
                            let id = db.insert_document(&imp.name, imp.kind, &imp.document, None, &desc, folder, &|f| {
                                p.set_fraction(base + (0.7 + 0.3 * f) / n)
                            })?;
                            ids.push(id);
                        }
                        Ok(ids)
                    },
                    |r| Message::ImportDone(Payload::new(r.unwrap_or_else(|| Err("Import failed".into())))),
                );
            }
            Message::ImportPicked(_) => {}
            Message::ImportDone(p) => {
                self.job = None;
                match p.take() {
                    Some(Ok(ids)) => {
                        self.reload_library();
                        self.status = format!("Imported {} file(s)", ids.len());
                        if let Some(&last) = ids.last() {
                            return self.open_document(last);
                        }
                    }
                    Some(Err(e)) => self.error = Some(e),
                    None => {}
                }
            }
            Message::AskRename(id) => {
                if let Some(d) = self.library.iter().find(|d| d.id == id) {
                    self.modal = Some(Modal::Rename { id, text: d.name.clone(), folder: false });
                    return dialogs::focus_rename();
                }
            }
            Message::AskDelete(id) => {
                if let Some(d) = self.library.iter().find(|d| d.id == id) {
                    self.modal = Some(Modal::ConfirmDelete { id, name: d.name.clone() });
                }
            }
            Message::ConfirmDelete(id) => {
                self.modal = None;
                if let Err(e) = self.db.delete_document(id) {
                    self.error = Some(e);
                }
                if self.open.as_ref().is_some_and(|o| o.info.id == id) {
                    self.open = None;
                }
                if self.library_selected == Some(id) {
                    self.library_selected = None;
                }
                self.reload_library();
            }
            Message::View(m) => self.on_view(m),
            Message::ZoomIn => self.zoom_by(1.5),
            Message::ZoomOut => self.zoom_by(1.0 / 1.5),
            Message::ZoomFit => self.with_doc(|od, _| {
                od.vp.col_w = od.vp.min_col_w(od.width());
                od.clamp_view();
            }),
            Message::ZoomLetters => self.with_doc(|od, _| {
                let center = od.vp.scroll_x + od.vp.visible_cols() / 2.0;
                od.vp.col_w = BASE_COL_W;
                od.vp.scroll_x = center - od.vp.visible_cols() / 2.0;
                od.clamp_view();
            }),
            Message::RowHeight(h) => {
                self.settings.row_height = h;
                self.with_doc(|od, _| od.set_row_h(h));
                self.save_settings();
            }
            Message::JumpInput(s) => self.with_doc(|od, _| od.jump_text = s),
            Message::JumpSubmit => self.jump(),
            Message::FilterInput(s) => self.with_doc(|od, _| {
                od.display.name_filter = s;
                od.rebuild_items();
            }),
            Message::SearchInput(s) => self.with_doc(|od, _| {
                od.search.query = s;
                od.search.searched = false;
            }),
            Message::SearchSubmit => {
                let searched = self.open.as_ref().is_some_and(|o| o.search.searched && !o.search.hits.is_empty());
                if searched {
                    self.search_step(1);
                } else {
                    self.run_search();
                }
            }
            Message::SearchStep(d) => self.search_step(d),
            Message::SearchScope(s) => self.with_doc(|od, _| {
                od.search.scope = s;
                od.search.searched = false;
            }),
            Message::SearchMismatches(m) => self.with_doc(|od, _| {
                od.search.mismatches = m;
                od.search.searched = false;
            }),
            Message::SetHighlight(h) => {
                self.set_pref(|p| p.highlight = h);
                let prefs = self.settings.view.clone();
                self.with_doc(|od, _| od.recompute_graphs(&prefs));
            }
            Message::ToggleDots(b) => self.set_pref(|p| p.use_dots = b),
            Message::ToggleGaps(b) => self.set_pref(|p| p.highlight_gaps = b),
            Message::ToggleConsensus(b) => self.set_pref(|p| p.show_consensus = b),
            Message::ToggleIdentity(b) => self.set_pref(|p| p.show_identity = b),
            Message::ToggleAnnotations(b) => self.set_pref(|p| p.show_annotations = b),
            Message::ToggleOverview(b) => self.set_pref(|p| p.show_overview = b),
            Message::ToggleColors(b) => self.set_pref(|p| p.color_bases = b),
            Message::SetGraph(g) => self.set_pref(|p| p.graph = g),
            Message::SetThreshold(t) => {
                self.set_pref(|p| p.threshold = t);
                let prefs = self.settings.view.clone();
                self.with_doc(|od, _| od.recompute_consensus(&prefs));
                self.after_selection_change();
            }
            Message::SortBy(k) => self.with_doc(|od, _| {
                od.display.sort = k;
                od.rebuild_items();
            }),
            Message::SortDescending(b) => self.with_doc(|od, _| {
                od.display.descending = b;
                od.rebuild_items();
            }),
            Message::GroupBy(g) => {
                self.with_doc(|od, _| {
                    od.display.group_by = match g {
                        GroupChoice::None => None,
                        GroupChoice::Key(k) => Some(k),
                    };
                    od.display.folded_groups.clear();
                    od.rebuild_items();
                });
                self.after_selection_change();
            }
            Message::Collapse(b) => self.with_doc(|od, _| {
                od.display.collapse_identical = b;
                od.rebuild_items();
            }),
            Message::SetReferenceFromSelection => {
                let row = self.open.as_ref().and_then(|od| {
                    if od.selected_rows.len() == 1 { od.selected_rows.iter().next().copied() } else { None }
                });
                if let Some(row) = row {
                    self.set_reference(Some(row));
                }
            }
            Message::ClearReference => self.set_reference(None),
            Message::ClearRowSelection => self.with_doc(|od, _| od.selected_rows.clear()),
            Message::ClearSelection => self.with_doc(|od, _| {
                od.selection = None;
                od.report = None;
                od.popup_open = false;
            }),
            Message::ShowPopup => self.with_doc(|od, _| od.popup_open = od.report.is_some()),
            Message::ClosePopup => self.with_doc(|od, _| od.popup_open = false),
            Message::PrimerSource(s) => {
                self.settings.primer.source = s;
                self.save_settings();
                self.after_selection_change();
            }
            Message::PrimerOrientation(o) => {
                self.settings.primer.orientation = o;
                self.save_settings();
                self.after_selection_change();
            }
            Message::CopyOligo => {
                if let Some(r) = self.open.as_ref().and_then(|o| o.report.as_ref()) {
                    self.status = "Oligo copied to the clipboard".into();
                    return iced::clipboard::write(String::from_utf8_lossy(&r.oligo).to_string());
                }
            }
            Message::CopyReport => {
                if let Some(od) = self.open.as_ref()
                    && let Some(r) = od.report.as_ref() {
                        let text = primer::report_text(r, od.ref_range(r.col_start, r.col_end));
                        self.status = "Report copied to the clipboard".into();
                        return iced::clipboard::write(text);
                    }
            }
            Message::AnnotationName(s) => self.with_doc(|od, _| od.annotation_name = s),
            Message::AnnotationKind(k) => self.with_doc(|od, _| od.annotation_kind = k),
            Message::AddAnnotation => self.add_annotation(),
            Message::DeleteAnnotation(id) => {
                if let Err(e) = self.db.delete_annotation(id) {
                    self.error = Some(e);
                }
                self.with_doc(|od, _| od.annotations.retain(|a| a.id != id));
            }
            Message::GotoAnnotation(i) => self.goto_annotation(i),
            Message::ImportMetadata => {
                if self.open.is_none() {
                    return Task::none();
                }
                return Task::perform(
                    rfd::AsyncFileDialog::new()
                        .set_title("Import metadata table")
                        .add_filter("Tables", &["csv", "tsv", "txt", "tab"])
                        .add_filter("All files", &["*"])
                        .pick_file(),
                    |f| Message::MetadataPicked(f.map(|f| f.path().to_path_buf())),
                );
            }
            Message::MetadataPicked(Some(path)) => self.import_metadata(&path),
            Message::MetadataPicked(None) => {}
            Message::OpenNameFields => {
                if self.open.is_some() {
                    self.modal = Some(Modal::NameFields { delimiter: "|".into(), keys: String::new() });
                }
            }
            Message::OpenExport => {
                if let Some(od) = &self.open {
                    self.modal = Some(Modal::Export(ExportDialog::new(od)));
                }
            }
            Message::Export(m) => {
                if let Some(Modal::Export(d)) = self.modal.as_mut()
                    && let Some(task) = d.update(m) {
                        return task;
                    }
            }
            Message::ExportPicked(Some(path)) => return self.export(path),
            Message::ExportPicked(None) => {}
            Message::ExportDone(p) => {
                self.job = None;
                match p.take() {
                    Some(Ok(s)) => self.status = s,
                    Some(Err(e)) => self.error = Some(e),
                    None => {}
                }
            }
            Message::OpenMapDialog => {
                if let Some(od) = &self.open {
                    self.modal = Some(Modal::Map(MapDialog::new(od, &self.settings.map, &self.library)));
                } else {
                    self.error = Some("Open the sequences to map first (select a document in the library).".into());
                }
            }
            Message::Map(m) => {
                if let Some(Modal::Map(d)) = self.modal.as_mut() {
                    d.update(m, self.open.as_ref());
                }
            }
            Message::RunMapping => return self.run_mapping(),
            Message::MapDone(p) => {
                self.job = None;
                match p.take() {
                    Some(Ok((id, report))) => {
                        self.reload_library();
                        self.status = report.clone();
                        self.modal = Some(Modal::Info { title: "Map to reference finished".into(), text: report });
                        return self.open_document(id);
                    }
                    Some(Err(e)) => self.error = Some(e),
                    None => {}
                }
            }
            Message::OpenSettings => self.modal = Some(Modal::Settings(SettingsDialog::new(&self.settings))),
            Message::Settings(m) => {
                if let Some(Modal::Settings(d)) = self.modal.as_mut() {
                    d.update(m);
                }
            }
            Message::SaveSettings => {
                if let Some(Modal::Settings(d)) = self.modal.take() {
                    d.apply(&mut self.settings);
                    self.save_settings();
                    self.after_selection_change();
                }
            }
            Message::Modal(m) => return self.on_modal(m),
            Message::CloseModal => self.modal = None,
            Message::CancelJob => {
                if let Some(j) = &self.job {
                    j.progress.cancel.store(true, std::sync::atomic::Ordering::Relaxed);
                }
            }
            Message::DemoData => return self.demo_data(),
            Message::ToggleMenu(id) => self.menu = if self.menu == Some(id) { None } else { Some(id) },
            Message::HoverMenu(id) => {
                if self.menu.is_some() {
                    self.menu = Some(id);
                }
            }
            Message::CloseMenu => self.menu = None,
            Message::MenuAction(m) => {
                self.menu = None;
                return self.update(*m);
            }
            Message::CloseDocument => {
                if self.open.take().is_some() {
                    self.status = "Document closed".into();
                }
            }
            Message::FocusJump => {
                if self.open.is_some() {
                    return Task::batch([iced::widget::operation::focus(JUMP_INPUT), iced::widget::operation::select_all(JUMP_INPUT)]);
                }
            }
            Message::FocusSearch => {
                if self.open.is_some() {
                    return Task::batch([iced::widget::operation::focus(SEARCH_INPUT), iced::widget::operation::select_all(SEARCH_INPUT)]);
                }
            }
            Message::Exit => return iced::exit(),
        }
        Task::none()
    }

    fn set_pref(&mut self, f: impl FnOnce(&mut ViewPrefs)) {
        f(&mut self.settings.view);
        self.save_settings();
    }

    fn zoom_by(&mut self, factor: f32) {
        self.with_doc(|od, _| {
            let center = od.vp.scroll_x + od.vp.visible_cols() / 2.0;
            od.vp.col_w *= factor;
            od.clamp_view();
            od.vp.scroll_x = center - od.vp.visible_cols() / 2.0;
            od.clamp_view();
        });
    }

    fn on_view(&mut self, m: ViewMsg) {
        let Some(od) = self.open.as_mut() else { return };
        let width = od.width();
        match m {
            ViewMsg::Resized { canvas, seq_w, rows_h } => {
                let first = od.vp.canvas == (0.0, 0.0);
                od.vp.canvas = canvas;
                od.vp.seq_w = seq_w;
                od.vp.rows_h = rows_h;
                if first && od.vp.col_w * width as f32 <= seq_w {
                    od.vp.col_w = od.vp.col_w.max(od.vp.min_col_w(width));
                }
                od.clamp_view();
            }
            ViewMsg::Scroll { dx_cols, dy_px } => {
                od.vp.scroll_x += dx_cols;
                od.vp.scroll_y += dy_px;
                od.clamp_view();
            }
            ViewMsg::SetScroll { x, y } => {
                if let Some(x) = x {
                    od.vp.scroll_x = x;
                }
                if let Some(y) = y {
                    od.vp.scroll_y = y;
                }
                od.clamp_view();
            }
            ViewMsg::Zoom { factor, anchor_col, anchor_frac } => {
                od.vp.col_w *= factor;
                od.clamp_view();
                od.vp.scroll_x = anchor_col - anchor_frac * od.vp.visible_cols();
                od.clamp_view();
            }
            ViewMsg::Select { anchor, to, finished } => {
                od.selection = Some(Selection::new(anchor, to));
                if !finished {
                    // The report belongs to the previous selection until the drag ends.
                    od.popup_open = false;
                }
                if finished {
                    let opts = self.settings.primer;
                    od.evaluate_selection(&opts);
                    od.popup_open = od.selection.is_some_and(|s| s.len() >= 1);
                    if od.annotation_name.is_empty() {
                        od.annotation_name = format!("Oligo {}", od.annotations.len() + 1);
                    }
                }
            }
            ViewMsg::ClickItem { item, ctrl, shift } => {
                let Some(row) = od.item_row(item) else { return };
                if shift {
                    if let Some(last) = od.last_clicked_item {
                        let (a, b) = (last.min(item), last.max(item));
                        for i in a..=b {
                            if let Some(r) = od.item_row(i) {
                                od.selected_rows.insert(r);
                            }
                        }
                    } else {
                        od.selected_rows.insert(row);
                    }
                } else if ctrl {
                    if !od.selected_rows.remove(&row) {
                        od.selected_rows.insert(row);
                    }
                } else if od.selected_rows.len() == 1 && od.selected_rows.contains(&row) {
                    od.selected_rows.clear();
                } else {
                    od.selected_rows.clear();
                    od.selected_rows.insert(row);
                }
                od.last_clicked_item = Some(item);
            }
            ViewMsg::ToggleGroup(label) => {
                if !od.display.folded_groups.remove(&label) {
                    od.display.folded_groups.insert(label);
                }
                od.rebuild_items();
            }
            ViewMsg::Hover(h) => od.hover = h,
            ViewMsg::ClickAnnotation(i) => self.goto_annotation(i),
        }
    }

    fn on_key(&mut self, e: keyboard::Event) -> Task<Message> {
        let keyboard::Event::KeyPressed { key, modifiers, .. } = e else { return Task::none() };
        if let Key::Named(Named::Escape) = key
            && self.menu.is_some()
        {
            self.menu = None;
            return Task::none();
        }
        if let Key::Named(Named::Escape) = key
            && self.drag.is_some()
        {
            self.drag = None;
            self.status = "Move cancelled".into();
            return Task::none();
        }
        if let Key::Named(Named::Escape) = key {
            if self.modal.is_some() {
                self.modal = None;
            } else if self.open.as_ref().is_some_and(|o| o.popup_open) {
                self.with_doc(|od, _| od.popup_open = false);
            } else {
                return self.update(Message::ClearSelection);
            }
            return Task::none();
        }
        if self.modal.is_some() {
            return Task::none();
        }
        if modifiers.command() {
            if let Key::Character(c) = &key {
                match c.as_str() {
                    "c" => return self.update(Message::CopyOligo),
                    "o" => return self.update(Message::ImportFasta),
                    "e" => return self.update(Message::OpenExport),
                    "m" => return self.update(Message::OpenMapDialog),
                    "w" => return self.update(Message::CloseDocument),
                    "g" => return self.update(Message::FocusJump),
                    "f" => return self.update(Message::FocusSearch),
                    "," => return self.update(Message::OpenSettings),
                    "=" | "+" => self.zoom_by(1.5),
                    "-" => self.zoom_by(1.0 / 1.5),
                    "0" => return self.update(Message::ZoomFit),
                    _ => {}
                }
            }
            return Task::none();
        }
        let Some(od) = self.open.as_mut() else { return Task::none() };
        let page_x = od.vp.visible_cols() * 0.9;
        let page_y = od.vp.rows_h * 0.9;
        let row_h = od.vp.row_h;
        let step_x = (1.0f32).max(24.0 / od.vp.col_w);
        let (dx, dy) = match key {
            Key::Named(Named::ArrowLeft) => (if modifiers.shift() { -page_x } else { -step_x }, 0.0),
            Key::Named(Named::ArrowRight) => (if modifiers.shift() { page_x } else { step_x }, 0.0),
            Key::Named(Named::ArrowUp) => (0.0, -row_h),
            Key::Named(Named::ArrowDown) => (0.0, row_h),
            Key::Named(Named::PageUp) => (0.0, -page_y),
            Key::Named(Named::PageDown) => (0.0, page_y),
            Key::Named(Named::Home) => (-(od.width() as f32), 0.0),
            Key::Named(Named::End) => (od.width() as f32, 0.0),
            _ => return Task::none(),
        };
        od.vp.scroll_x += dx;
        od.vp.scroll_y += dy;
        od.clamp_view();
        Task::none()
    }

    fn jump(&mut self) {
        let Some(od) = self.open.as_mut() else { return };
        let text: String = od.jump_text.chars().filter(|c| c.is_ascii_digit() || *c == '-').collect();
        let mut parts = text.split('-').filter(|s| !s.is_empty()).map(|s| s.parse::<usize>().ok());
        let Some(Some(a)) = parts.next() else {
            self.error = Some("Enter a position like 19160 or a range like 19160-19180".into());
            return;
        };
        let b = parts.next().flatten();
        let to_col = |pos: usize| -> usize {
            match &od.ref_coords {
                Some(rc) => rc.column_of(pos).unwrap_or(0),
                None => pos.clamp(1, od.doc.width.max(1)) - 1,
            }
        };
        let c0 = to_col(a);
        let c1 = b.map(to_col).unwrap_or(c0);
        let (c0, c1) = (c0.min(c1), c0.max(c1));
        if b.is_some() {
            od.selection = Some(Selection::new(c0, c1));
            let opts = self.settings.primer;
            od.evaluate_selection(&opts);
            od.popup_open = true;
        }
        if od.vp.col_w < 4.0 {
            od.vp.col_w = BASE_COL_W;
        }
        od.clamp_view();
        od.vp.scroll_x = (c0 + c1) as f32 / 2.0 - od.vp.visible_cols() / 2.0;
        od.clamp_view();
    }

    fn run_search(&mut self) {
        let Some(od) = self.open.as_mut() else { return };
        let q = od.search.query.trim().to_string();
        if q.is_empty() {
            od.search.hits.clear();
            od.search.current = None;
            return;
        }
        let hits = search::search(&od.doc, od.reference, &od.consensus, &q, od.search.scope, od.search.mismatches as usize, true, 100_000);
        od.search.searched = true;
        od.search.current = if hits.is_empty() { None } else { Some(0) };
        self.status = match hits.len() {
            0 => format!("No match for {q}"),
            100_000 => "100,000+ matches (showing the first 100,000)".into(),
            n => format!("{} match(es) for {q}", canvas::group_digits(n)),
        };
        od.search.hits = hits;
        self.search_step(0);
    }

    fn search_step(&mut self, delta: i32) {
        let opts = self.settings.primer;
        let Some(od) = self.open.as_mut() else { return };
        if od.search.hits.is_empty() {
            return;
        }
        let n = od.search.hits.len() as i32;
        let cur = od.search.current.map(|c| c as i32).unwrap_or(0);
        let next = (cur + delta).rem_euclid(n) as usize;
        od.search.current = Some(next);
        let h = od.search.hits[next];
        if od.vp.col_w < 4.0 {
            od.vp.col_w = BASE_COL_W;
            od.clamp_view();
        }
        od.vp.reveal_col(h.c0, h.c1);
        if let Some(item) = h.row.and_then(|r| od.item_of_row(r)) {
            od.reveal_item(item);
        }
        od.clamp_view();
        od.selection = Some(Selection { anchor: h.c0, c0: h.c0, c1: h.c1 });
        od.evaluate_selection(&opts);
    }

    fn set_reference(&mut self, reference: Option<usize>) {
        let prefs = self.settings.view.clone();
        let opts = self.settings.primer;
        let Some(od) = self.open.as_mut() else { return };
        od.reference = reference;
        od.info.reference = reference;
        if let Err(e) = self.db.set_reference(od.info.id, reference) {
            self.error = Some(e);
        }
        od.recompute_stats(&prefs);
        od.selected_rows.clear();
        od.rebuild_items();
        od.evaluate_selection(&opts);
        self.reload_library();
        self.status = match reference {
            Some(r) => format!("Reference set to {}", self.open.as_ref().unwrap().doc.rows[r].name),
            None => "Reference cleared".into(),
        };
    }

    fn add_annotation(&mut self) {
        let Some(od) = self.open.as_mut() else { return };
        let Some(sel) = od.selection else { return };
        let name = if od.annotation_name.trim().is_empty() { format!("Oligo {}", od.annotations.len() + 1) } else { od.annotation_name.trim().to_string() };
        let note = od.report.as_ref().map(|r| String::from_utf8_lossy(&r.oligo).to_string()).unwrap_or_default();
        let mut a = Annotation { id: 0, name, kind: od.annotation_kind, start: sel.c0, end: sel.c1, note };
        match self.db.add_annotation(od.info.id, &a) {
            Ok(id) => {
                a.id = id;
                self.status = format!("Annotation “{}” added", a.name);
                od.annotations.push(a);
                od.annotations.sort_by_key(|a| (a.start, a.end));
                od.annotation_name = format!("Oligo {}", od.annotations.len() + 1);
            }
            Err(e) => self.error = Some(e),
        }
    }

    fn goto_annotation(&mut self, i: usize) {
        let opts = self.settings.primer;
        let Some(od) = self.open.as_mut() else { return };
        let Some(a) = od.annotations.get(i).cloned() else { return };
        od.selection = Some(Selection { anchor: a.start, c0: a.start, c1: a.end });
        match a.kind {
            AnnotationKind::ReversePrimer => self.settings.primer.orientation = Orientation::Reverse,
            AnnotationKind::ForwardPrimer => self.settings.primer.orientation = Orientation::Forward,
            _ => {}
        }
        let opts = PrimerOptions { orientation: self.settings.primer.orientation, ..opts };
        od.evaluate_selection(&opts);
        od.popup_open = true;
        od.annotation_name = a.name.clone();
        od.annotation_kind = a.kind;
        if od.vp.col_w < 4.0 {
            od.vp.col_w = BASE_COL_W;
            od.clamp_view();
        }
        od.vp.reveal_col(a.start, a.end);
        od.clamp_view();
    }

    fn apply_metadata(&mut self, f: impl FnOnce(&mut [Row]) -> String) {
        let Some(od) = self.open.as_mut() else { return };
        let doc = Arc::make_mut(&mut od.doc);
        let msg = f(&mut doc.rows);
        if let Err(e) = self.db.update_metadata(od.info.id, &doc.rows) {
            self.error = Some(e);
            return;
        }
        od.meta_keys = od.doc.meta_keys();
        od.rebuild_items();
        self.status = msg;
    }

    fn import_metadata(&mut self, path: &std::path::Path) {
        let table = match fasta::read_meta_table(path) {
            Ok(t) => t,
            Err(e) => {
                self.error = Some(e);
                return;
            }
        };
        self.apply_metadata(|rows| {
            let mut by_name: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
            for (i, (name, _)) in table.rows.iter().enumerate() {
                by_name.insert(name.to_lowercase(), i);
            }
            let mut matched = 0;
            for r in rows.iter_mut() {
                let key = r.name.to_lowercase();
                let alt = key.split('.').next().unwrap_or("").to_string();
                let idx = by_name.get(&key).or_else(|| by_name.get(&alt)).copied();
                if let Some(i) = idx {
                    matched += 1;
                    for (k, v) in table.keys.iter().zip(table.rows[i].1.iter()) {
                        if !v.is_empty() {
                            r.set_meta(k, v.clone());
                        }
                    }
                }
            }
            format!("Metadata: {} of {} sequences matched ({} columns)", matched, rows.len(), table.keys.len())
        });
    }

    fn on_modal(&mut self, m: dialogs::ModalMsg) -> Task<Message> {
        use dialogs::ModalMsg;
        match m {
            ModalMsg::RenameInput(s) => {
                if let Some(Modal::Rename { text, .. }) = self.modal.as_mut() {
                    *text = s;
                }
            }
            ModalMsg::RenameConfirm => {
                if let Some(Modal::Rename { id, text, folder }) = self.modal.take() {
                    let name = text.trim();
                    if !name.is_empty() {
                        let r = if folder { self.db.rename_folder(id, name) } else { self.db.rename_document(id, name) };
                        if let Err(e) = r {
                            self.error = Some(e);
                        }
                        if let Some(od) = self.open.as_mut().filter(|o| !folder && o.info.id == id) {
                            od.info.name = name.to_string();
                        }
                        self.reload_library();
                    }
                }
            }
            ModalMsg::DelimiterInput(s) => {
                if let Some(Modal::NameFields { delimiter, .. }) = self.modal.as_mut() {
                    *delimiter = s;
                }
            }
            ModalMsg::KeysInput(s) => {
                if let Some(Modal::NameFields { keys, .. }) = self.modal.as_mut() {
                    *keys = s;
                }
            }
            ModalMsg::ApplyNameFields => {
                if let Some(Modal::NameFields { delimiter, keys }) = self.modal.take() {
                    let delim = delimiter.chars().next().unwrap_or('|');
                    let names: Vec<String> = keys.split(',').map(|k| k.trim().to_string()).collect();
                    self.apply_metadata(|rows| {
                        let mut max_fields = 0;
                        for r in rows.iter_mut() {
                            let fields = fasta::name_fields(&r.name, delim);
                            if fields.len() < 2 {
                                continue;
                            }
                            max_fields = max_fields.max(fields.len());
                            for (i, v) in fields.into_iter().enumerate() {
                                let key = names.get(i).filter(|k| !k.is_empty()).cloned().unwrap_or_else(|| format!("Name field {}", i + 1));
                                r.set_meta(&key, v);
                            }
                        }
                        format!("Split sequence names into {max_fields} metadata field(s)")
                    });
                }
            }
            ModalMsg::DismissError => self.error = None,
        }
        Task::none()
    }

    fn export(&mut self, path: PathBuf) -> Task<Message> {
        let Some(Modal::Export(d)) = self.modal.take() else { return Task::none() };
        let Some(od) = self.open.as_ref() else { return Task::none() };
        let doc = od.doc.clone();
        let reference = od.reference;
        let consensus = od.consensus.clone();
        let region = match (d.selection_only, od.selection) {
            (true, Some(s)) => (s.c0, s.c1),
            _ => (0, doc.width),
        };
        let rows: Vec<usize> = match d.rows {
            ExportRows::All => (0..doc.rows.len()).filter(|&i| Some(i) != reference).collect(),
            ExportRows::Selected => {
                let mut v: Vec<usize> = od.selected_rows.iter().copied().collect();
                v.sort();
                v
            }
            ExportRows::Shown => od
                .items
                .iter()
                .filter_map(|it| match it {
                    crate::display::Item::Seq { row, .. } => Some(*row),
                    _ => None,
                })
                .collect(),
        };
        self.start_job(
            "Exporting FASTA",
            false,
            move |_p| {
                let (c0, c1) = region;
                let mut records = Vec::new();
                let render = |data: Vec<u8>| -> Vec<u8> {
                    let v: Vec<u8> = data.into_iter().map(|b| if b == b' ' { crate::seq::GAP } else { b }).collect();
                    if d.gapped { v } else { crate::seq::ungap(&v) }
                };
                if d.include_consensus {
                    records.push(fasta::Record { name: "Consensus".into(), description: String::new(), seq: render(consensus[c0..c1].to_vec()) });
                }
                if let (true, Some(r)) = (d.include_reference, reference) {
                    let row = &doc.rows[r];
                    records.push(fasta::Record { name: row.name.clone(), description: row.description.clone(), seq: render(row.slice_with_blanks(c0, c1)) });
                }
                for i in rows {
                    let row = &doc.rows[i];
                    let seq = render(row.slice_with_blanks(c0, c1));
                    if seq.iter().all(|&b| b == crate::seq::GAP) {
                        continue;
                    }
                    records.push(fasta::Record { name: row.name.clone(), description: row.description.clone(), seq });
                }
                fasta::write_fasta(&path, &records)?;
                Ok(format!("Exported {} sequences to {}", records.len(), path.display()))
            },
            |r| Message::ExportDone(Payload::new(r.unwrap_or_else(|| Err("Export failed".into())))),
        )
    }

    fn run_mapping(&mut self) -> Task<Message> {
        let Some(Modal::Map(d)) = self.modal.as_ref() else { return Task::none() };
        let Some(od) = self.open.as_ref() else { return Task::none() };
        let params = match d.params() {
            Ok(p) => p,
            Err(e) => {
                self.error = Some(e);
                return Task::none();
            }
        };
        let choice = d.reference.clone();
        let name = if d.name.trim().is_empty() { format!("{} mapped", od.info.name) } else { d.name.trim().to_string() };
        let only_selected = d.only_selected && !od.selected_rows.is_empty();
        let save_unmapped = d.save_unmapped;
        let doc = od.doc.clone();
        let ref_row_this = match &choice {
            RefChoice::ThisDocument(r) => *r,
            RefChoice::OtherDocument(_) => None,
        };
        let mut query_rows: Vec<usize> = if only_selected {
            let mut v: Vec<usize> = od.selected_rows.iter().copied().collect();
            v.sort();
            v
        } else {
            (0..doc.rows.len()).collect()
        };
        query_rows.retain(|&i| Some(i) != ref_row_this);
        if od.info.kind == DocKind::Contig
            && let Some(r) = od.reference {
                query_rows.retain(|&i| i != r);
            }
        let od_folder = od.info.folder;
        self.settings.map = params.clone();
        self.save_settings();
        self.modal = None;
        let db_path = self.db_path.clone();
        let folder = od_folder;
        self.start_job(
            "Map to reference",
            true,
            move |p| {
                let mut db = Db::open(&db_path)?;
                let reference: Row = match choice {
                    RefChoice::ThisDocument(Some(r)) => doc.rows[r].clone(),
                    RefChoice::ThisDocument(None) => return Err("Choose a reference sequence".to_string()),
                    RefChoice::OtherDocument(id) => {
                        p.set(0.0, "Loading reference");
                        let (info, rdoc, _) = db.load_document(id, &|_| {})?;
                        let r = info.reference.unwrap_or(0);
                        rdoc.rows.get(r).cloned().ok_or("The reference document is empty")?
                    }
                };
                let queries: Vec<&Row> = query_rows.iter().map(|&i| &doc.rows[i]).collect();
                if queries.is_empty() {
                    return Err("There are no sequences to map".into());
                }
                let progress = |f: f32, m: &str| p.set(f * 0.95, m);
                let ctl = mapper::Control { progress: &progress, cancel: &p.cancel };
                let out = mapper::map_to_reference(&reference, &queries, &params, &ctl)?;
                let r = &out.report;
                let report = format!(
                    "{} of {} sequences mapped to {} ({} reverse complemented), mean identity {:.2}%, {} iteration(s), {:.1} s",
                    canvas::group_digits(r.mapped),
                    canvas::group_digits(r.queries),
                    reference.name,
                    canvas::group_digits(r.reverse),
                    r.mean_identity,
                    r.iterations,
                    r.seconds
                );
                p.set(0.96, "Saving contig");
                let id = db.insert_document(&name, DocKind::Contig, &out.contig, Some(0), &report, folder, &|f| p.set_fraction(0.96 + f * 0.03))?;
                if save_unmapped && !out.unmapped.is_empty() {
                    let un = Document::new(out.unmapped);
                    db.insert_document(&format!("{name} – unmapped"), DocKind::Sequences, &un, None, "Sequences that did not map", folder, &|_| {})?;
                }
                Ok((id, report))
            },
            |r| Message::MapDone(Payload::new(r.unwrap_or_else(|| Err("Mapping failed unexpectedly".into())))),
        )
    }

    /// Saves the selected columns of all sequences (without the reference) as a
    /// new library document next to the open one.
    fn extract_slice(&mut self) -> Task<Message> {
        if self.busy() {
            return Task::none();
        }
        let Some(od) = self.open.as_ref() else { return Task::none() };
        let Some(sel) = od.selection else {
            self.error = Some("Select a column range first (drag across the alignment).".into());
            return Task::none();
        };
        let doc = od.doc.clone();
        let reference = od.reference;
        let folder = od.info.folder;
        let kind = if od.info.kind == DocKind::Sequences { DocKind::Sequences } else { DocKind::Alignment };
        let range = match od.ref_range(sel.c0, sel.c1) {
            Some((a, b)) => format!("ref {a}-{b}"),
            None => format!("columns {}-{}", sel.c0 + 1, sel.c1),
        };
        let name = format!("{} [{range}]", od.info.name);
        let source = od.info.name.clone();
        let db_path = self.db_path.clone();
        self.start_job(
            "Extracting slice",
            false,
            move |p| {
                let slice = doc.slice_columns(sel.c0, sel.c1, reference);
                if slice.rows.is_empty() {
                    return Err("No sequence has bases in the selected columns".to_string());
                }
                let mut db = Db::open(&db_path)?;
                let desc = format!("Slice {range} of “{source}”");
                let id = db.insert_document(&name, kind, &slice, None, &desc, folder, &|f| p.set_fraction(f))?;
                Ok((id, format!("Saved “{name}”: {} sequences × {} columns", canvas::group_digits(slice.rows.len()), slice.width)))
            },
            |r| Message::SliceDone(Payload::new(r.unwrap_or_else(|| Err("Extracting the slice failed".into())))),
        )
    }

    fn demo_data(&mut self) -> Task<Message> {
        if self.busy() {
            return Task::none();
        }
        let db_path = self.db_path.clone();
        let folder = self.target_folder();
        self.start_job(
            "Generating demo data",
            false,
            move |p| {
                let (reference, queries) = crate::demo::generate(35_000, 2_000, 7);
                let mut db = Db::open(&db_path)?;
                p.set(0.3, "Saving demo sequences");
                let mut rows = vec![reference];
                rows.extend(queries);
                let doc = Document::new(rows);
                let id = db.insert_document("Demo genomes (35 kb × 2,000)", DocKind::Sequences, &doc, Some(0), "Synthetic demo data", folder, &|f| {
                    p.set_fraction(0.3 + 0.7 * f)
                })?;
                Ok(vec![id])
            },
            |r| Message::ImportDone(Payload::new(r.unwrap_or_else(|| Err("Demo generation failed".into())))),
        )
    }
}
