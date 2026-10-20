//! SQLite document library, stored next to the executable (`pcr_studio.db`).

use crate::model::{Annotation, AnnotationKind, DocInfo, DocKind, Document, Folder, Row};
use rusqlite::{Connection, OptionalExtension, params};
use std::path::{Path, PathBuf};

pub const DB_FILE_NAME: &str = "pcr_studio.db";

/// Location of the database: the directory of the running executable.
pub fn default_path() -> PathBuf {
    std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|d| d.to_path_buf()))
        .unwrap_or_else(|| PathBuf::from("."))
        .join(DB_FILE_NAME)
}

pub struct Db {
    conn: Connection,
}

const SCHEMA: &str = "
CREATE TABLE IF NOT EXISTS documents (
    id          INTEGER PRIMARY KEY,
    name        TEXT NOT NULL,
    kind        TEXT NOT NULL,
    n_rows      INTEGER NOT NULL,
    width       INTEGER NOT NULL,
    reference   INTEGER,
    description TEXT NOT NULL DEFAULT '',
    created     TEXT NOT NULL DEFAULT (datetime('now', 'localtime'))
);
CREATE TABLE IF NOT EXISTS rows (
    doc_id      INTEGER NOT NULL REFERENCES documents(id) ON DELETE CASCADE,
    idx         INTEGER NOT NULL,
    name        TEXT NOT NULL,
    description TEXT NOT NULL,
    start       INTEGER NOT NULL,
    data        BLOB NOT NULL,
    meta        TEXT NOT NULL,
    PRIMARY KEY (doc_id, idx)
);
CREATE TABLE IF NOT EXISTS annotations (
    id          INTEGER PRIMARY KEY,
    doc_id      INTEGER NOT NULL REFERENCES documents(id) ON DELETE CASCADE,
    name        TEXT NOT NULL,
    kind        TEXT NOT NULL,
    start_col   INTEGER NOT NULL,
    end_col     INTEGER NOT NULL,
    note        TEXT NOT NULL DEFAULT ''
);
CREATE TABLE IF NOT EXISTS folders (
    id          INTEGER PRIMARY KEY,
    name        TEXT NOT NULL,
    parent_id   INTEGER REFERENCES folders(id) ON DELETE CASCADE,
    collapsed   INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS settings (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS texts (
    doc_id      INTEGER PRIMARY KEY REFERENCES documents(id) ON DELETE CASCADE,
    content     TEXT NOT NULL
);
";

fn meta_to_json(meta: &[(String, String)]) -> String {
    serde_json::to_string(meta).unwrap_or_else(|_| "[]".into())
}

fn meta_from_json(s: &str) -> Vec<(String, String)> {
    serde_json::from_str(s).unwrap_or_default()
}

impl Db {
    pub fn open(path: &Path) -> Result<Db, String> {
        let conn = Connection::open(path).map_err(|e| format!("Cannot open database {}: {e}", path.display()))?;
        conn.execute_batch("PRAGMA foreign_keys = ON; PRAGMA journal_mode = WAL; PRAGMA synchronous = NORMAL;")
            .map_err(|e| e.to_string())?;
        conn.execute_batch(SCHEMA).map_err(|e| format!("Cannot initialize database: {e}"))?;
        // Libraries created before folders existed get the column added.
        let has_folder: bool = conn
            .prepare("SELECT 1 FROM pragma_table_info('documents') WHERE name = 'folder_id'")
            .and_then(|mut st| st.exists([]))
            .map_err(|e| e.to_string())?;
        if !has_folder {
            conn.execute_batch("ALTER TABLE documents ADD COLUMN folder_id INTEGER REFERENCES folders(id) ON DELETE SET NULL;")
                .map_err(|e| format!("Cannot upgrade database: {e}"))?;
        }
        Ok(Db { conn })
    }

    pub fn list_documents(&self) -> Result<Vec<DocInfo>, String> {
        let mut st = self
            .conn
            .prepare("SELECT id, name, kind, n_rows, width, reference, created, description, folder_id FROM documents ORDER BY id")
            .map_err(|e| e.to_string())?;
        let rows = st
            .query_map([], |r| {
                Ok(DocInfo {
                    id: r.get(0)?,
                    name: r.get(1)?,
                    kind: DocKind::parse(&r.get::<_, String>(2)?),
                    n_rows: r.get::<_, i64>(3)? as usize,
                    width: r.get::<_, i64>(4)? as usize,
                    reference: r.get::<_, Option<i64>>(5)?.map(|x| x as usize),
                    created: r.get(6)?,
                    description: r.get(7)?,
                    folder: r.get(8)?,
                })
            })
            .map_err(|e| e.to_string())?;
        rows.collect::<Result<Vec<_>, _>>().map_err(|e| e.to_string())
    }

    pub fn document_info(&self, id: i64) -> Result<Option<DocInfo>, String> {
        Ok(self.list_documents()?.into_iter().find(|d| d.id == id))
    }

    /// Stores a document; returns its id.
    pub fn insert_document(
        &mut self,
        name: &str,
        kind: DocKind,
        doc: &Document,
        reference: Option<usize>,
        description: &str,
        folder: Option<i64>,
        progress: &dyn Fn(f32),
    ) -> Result<i64, String> {
        let tx = self.conn.transaction().map_err(|e| e.to_string())?;
        tx.execute(
            "INSERT INTO documents (name, kind, n_rows, width, reference, description, folder_id) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
            params![name, kind.as_str(), doc.rows.len() as i64, doc.width as i64, reference.map(|r| r as i64), description, folder],
        )
        .map_err(|e| e.to_string())?;
        let id = tx.last_insert_rowid();
        {
            let mut st = tx
                .prepare("INSERT INTO rows (doc_id, idx, name, description, start, data, meta) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)")
                .map_err(|e| e.to_string())?;
            let n = doc.rows.len().max(1);
            for (i, r) in doc.rows.iter().enumerate() {
                st.execute(params![id, i as i64, r.name, r.description, r.start as i64, r.data, meta_to_json(&r.meta)])
                    .map_err(|e| e.to_string())?;
                if i % 500 == 0 {
                    progress(i as f32 / n as f32);
                }
            }
        }
        tx.commit().map_err(|e| e.to_string())?;
        Ok(id)
    }

    /// Stores a text document (no sequences); returns its id. A folder that no
    /// longer exists is replaced by the top level.
    pub fn insert_text_document(&mut self, name: &str, content: &str, description: &str, folder: Option<i64>) -> Result<i64, String> {
        let tx = self.conn.transaction().map_err(|e| e.to_string())?;
        let folder = match folder {
            Some(f) => tx.query_row("SELECT id FROM folders WHERE id = ?1", params![f], |r| r.get::<_, i64>(0)).optional().map_err(|e| e.to_string())?,
            None => None,
        };
        tx.execute(
            "INSERT INTO documents (name, kind, n_rows, width, reference, description, folder_id) VALUES (?1, ?2, 0, 0, NULL, ?3, ?4)",
            params![name, DocKind::Text.as_str(), description, folder],
        )
        .map_err(|e| e.to_string())?;
        let id = tx.last_insert_rowid();
        tx.execute("INSERT INTO texts (doc_id, content) VALUES (?1, ?2)", params![id, content]).map_err(|e| e.to_string())?;
        tx.commit().map_err(|e| e.to_string())?;
        Ok(id)
    }

    /// Content of a text document.
    pub fn text_content(&self, id: i64) -> Result<String, String> {
        self.conn
            .query_row("SELECT content FROM texts WHERE doc_id = ?1", params![id], |r| r.get(0))
            .optional()
            .map_err(|e| e.to_string())?
            .ok_or_else(|| format!("Text document {id} not found"))
    }

    /// Local date and time as the library records it (`YYYY-MM-DD HH:MM:SS`).
    pub fn now_local(&self) -> String {
        self.conn.query_row("SELECT datetime('now', 'localtime')", [], |r| r.get(0)).unwrap_or_default()
    }

    pub fn load_document(&self, id: i64, progress: &dyn Fn(f32)) -> Result<(DocInfo, Document, Vec<Annotation>), String> {
        let info = self.document_info(id)?.ok_or_else(|| format!("Document {id} not found"))?;
        let mut st = self
            .conn
            .prepare("SELECT name, description, start, data, meta FROM rows WHERE doc_id = ?1 ORDER BY idx")
            .map_err(|e| e.to_string())?;
        let mut q = st.query(params![id]).map_err(|e| e.to_string())?;
        let mut rows = Vec::with_capacity(info.n_rows);
        let n = info.n_rows.max(1);
        while let Some(r) = q.next().map_err(|e| e.to_string())? {
            rows.push(Row {
                name: r.get(0).map_err(|e| e.to_string())?,
                description: r.get(1).map_err(|e| e.to_string())?,
                start: r.get::<_, i64>(2).map_err(|e| e.to_string())? as usize,
                data: r.get(3).map_err(|e| e.to_string())?,
                meta: meta_from_json(&r.get::<_, String>(4).map_err(|e| e.to_string())?),
            });
            if rows.len() % 500 == 0 {
                progress(rows.len() as f32 / n as f32);
            }
        }
        let annotations = self.annotations(id)?;
        let mut doc = Document::new(rows);
        doc.width = doc.width.max(info.width);
        Ok((info, doc, annotations))
    }

    pub fn delete_document(&self, id: i64) -> Result<(), String> {
        self.conn.execute("DELETE FROM documents WHERE id = ?1", params![id]).map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn rename_document(&self, id: i64, name: &str) -> Result<(), String> {
        self.conn.execute("UPDATE documents SET name = ?2 WHERE id = ?1", params![id, name]).map_err(|e| e.to_string())?;
        Ok(())
    }

    /// Replaces the sequence of a single-sequence document and renames the
    /// document. Annotations that end past the new sequence are removed;
    /// returns how many.
    pub fn replace_single_row(&mut self, id: i64, doc_name: &str, kind: DocKind, row: &Row) -> Result<usize, String> {
        let tx = self.conn.transaction().map_err(|e| e.to_string())?;
        let n_rows: i64 = tx
            .query_row("SELECT n_rows FROM documents WHERE id = ?1", params![id], |r| r.get(0))
            .optional()
            .map_err(|e| e.to_string())?
            .ok_or_else(|| format!("Document {id} not found"))?;
        if n_rows != 1 {
            return Err("Only documents with a single sequence can be edited".into());
        }
        let width = row.end() as i64;
        tx.execute("UPDATE documents SET name = ?2, kind = ?3, width = ?4 WHERE id = ?1", params![id, doc_name, kind.as_str(), width])
            .map_err(|e| e.to_string())?;
        tx.execute(
            "UPDATE rows SET name = ?2, description = ?3, start = ?4, data = ?5, meta = ?6 WHERE doc_id = ?1 AND idx = 0",
            params![id, row.name, row.description, row.start as i64, row.data, meta_to_json(&row.meta)],
        )
        .map_err(|e| e.to_string())?;
        let dropped = tx.execute("DELETE FROM annotations WHERE doc_id = ?1 AND end_col > ?2", params![id, width]).map_err(|e| e.to_string())?;
        tx.commit().map_err(|e| e.to_string())?;
        Ok(dropped)
    }

    /// Stores a new row order: row `order[i]` becomes row `i`; `reference` is
    /// the reference's new index. Rows move cycle by cycle through one free
    /// slot, so each row is rewritten once (plus one extra move per cycle).
    pub fn reorder_rows(&mut self, id: i64, order: &[usize], reference: Option<usize>) -> Result<(), String> {
        let tx = self.conn.transaction().map_err(|e| e.to_string())?;
        let n_rows: i64 = tx
            .query_row("SELECT n_rows FROM documents WHERE id = ?1", params![id], |r| r.get(0))
            .optional()
            .map_err(|e| e.to_string())?
            .ok_or_else(|| format!("Document {id} not found"))?;
        let mut seen = vec![false; order.len()];
        if n_rows as usize != order.len() || order.iter().any(|&o| o >= order.len() || std::mem::replace(&mut seen[o], true)) {
            return Err("The new row order does not match the document".into());
        }
        {
            let mut mv = tx.prepare("UPDATE rows SET idx = ?3 WHERE doc_id = ?1 AND idx = ?2").map_err(|e| e.to_string())?;
            let mut mv = |from: i64, to: i64| mv.execute(params![id, from, to]).map(|_| ()).map_err(|e| e.to_string());
            const FREE: i64 = -1;
            let mut done = vec![false; order.len()];
            for s in 0..order.len() {
                if done[s] || order[s] == s {
                    continue;
                }
                // Park the row of slot `s`, then fill each emptied slot with
                // the row that belongs there until the cycle closes.
                mv(s as i64, FREE)?;
                let mut slot = s;
                loop {
                    done[slot] = true;
                    let from = order[slot];
                    if from == s {
                        mv(FREE, slot as i64)?;
                        break;
                    }
                    mv(from as i64, slot as i64)?;
                    slot = from;
                }
            }
        }
        tx.execute("UPDATE documents SET reference = ?2 WHERE id = ?1", params![id, reference.map(|r| r as i64)])
            .map_err(|e| e.to_string())?;
        tx.commit().map_err(|e| e.to_string())
    }

    pub fn set_reference(&self, id: i64, reference: Option<usize>) -> Result<(), String> {
        self.conn
            .execute("UPDATE documents SET reference = ?2 WHERE id = ?1", params![id, reference.map(|r| r as i64)])
            .map_err(|e| e.to_string())?;
        Ok(())
    }

    /// Replaces the metadata of all rows (after importing a table / parsing names).
    pub fn update_metadata(&mut self, id: i64, rows: &[Row]) -> Result<(), String> {
        let tx = self.conn.transaction().map_err(|e| e.to_string())?;
        {
            let mut st = tx.prepare("UPDATE rows SET meta = ?3 WHERE doc_id = ?1 AND idx = ?2").map_err(|e| e.to_string())?;
            for (i, r) in rows.iter().enumerate() {
                st.execute(params![id, i as i64, meta_to_json(&r.meta)]).map_err(|e| e.to_string())?;
            }
        }
        tx.commit().map_err(|e| e.to_string())
    }

    pub fn annotations(&self, doc_id: i64) -> Result<Vec<Annotation>, String> {
        let mut st = self
            .conn
            .prepare("SELECT id, name, kind, start_col, end_col, note FROM annotations WHERE doc_id = ?1 ORDER BY start_col, id")
            .map_err(|e| e.to_string())?;
        let rows = st
            .query_map(params![doc_id], |r| {
                Ok(Annotation {
                    id: r.get(0)?,
                    name: r.get(1)?,
                    kind: AnnotationKind::parse(&r.get::<_, String>(2)?),
                    start: r.get::<_, i64>(3)? as usize,
                    end: r.get::<_, i64>(4)? as usize,
                    note: r.get(5)?,
                })
            })
            .map_err(|e| e.to_string())?;
        rows.collect::<Result<Vec<_>, _>>().map_err(|e| e.to_string())
    }

    pub fn add_annotation(&self, doc_id: i64, a: &Annotation) -> Result<i64, String> {
        self.conn
            .execute(
                "INSERT INTO annotations (doc_id, name, kind, start_col, end_col, note) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
                params![doc_id, a.name, a.kind.as_str(), a.start as i64, a.end as i64, a.note],
            )
            .map_err(|e| e.to_string())?;
        Ok(self.conn.last_insert_rowid())
    }

    /// Saves name, type, range and note of an existing annotation.
    pub fn update_annotation(&self, a: &Annotation) -> Result<(), String> {
        self.conn
            .execute(
                "UPDATE annotations SET name = ?2, kind = ?3, start_col = ?4, end_col = ?5, note = ?6 WHERE id = ?1",
                params![a.id, a.name, a.kind.as_str(), a.start as i64, a.end as i64, a.note],
            )
            .map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn delete_annotation(&self, id: i64) -> Result<(), String> {
        self.conn.execute("DELETE FROM annotations WHERE id = ?1", params![id]).map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn list_folders(&self) -> Result<Vec<Folder>, String> {
        let mut st = self
            .conn
            .prepare("SELECT id, name, parent_id, collapsed FROM folders ORDER BY name COLLATE NOCASE, id")
            .map_err(|e| e.to_string())?;
        let rows = st
            .query_map([], |r| Ok(Folder { id: r.get(0)?, name: r.get(1)?, parent: r.get(2)?, collapsed: r.get::<_, i64>(3)? != 0 }))
            .map_err(|e| e.to_string())?;
        rows.collect::<Result<Vec<_>, _>>().map_err(|e| e.to_string())
    }

    pub fn create_folder(&self, name: &str, parent: Option<i64>) -> Result<i64, String> {
        self.conn
            .execute("INSERT INTO folders (name, parent_id) VALUES (?1, ?2)", params![name, parent])
            .map_err(|e| e.to_string())?;
        Ok(self.conn.last_insert_rowid())
    }

    pub fn rename_folder(&self, id: i64, name: &str) -> Result<(), String> {
        self.conn.execute("UPDATE folders SET name = ?2 WHERE id = ?1", params![id, name]).map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn set_folder_collapsed(&self, id: i64, collapsed: bool) -> Result<(), String> {
        self.conn
            .execute("UPDATE folders SET collapsed = ?2 WHERE id = ?1", params![id, collapsed as i64])
            .map_err(|e| e.to_string())?;
        Ok(())
    }

    /// Moves a folder under `parent` (`None` = top level). Refuses to create cycles.
    pub fn move_folder(&self, id: i64, parent: Option<i64>) -> Result<(), String> {
        let folders = self.list_folders()?;
        let mut cur = parent;
        while let Some(p) = cur {
            if p == id {
                return Err("A folder cannot be moved into itself".into());
            }
            cur = folders.iter().find(|f| f.id == p).and_then(|f| f.parent);
        }
        self.conn.execute("UPDATE folders SET parent_id = ?2 WHERE id = ?1", params![id, parent]).map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn move_document(&self, id: i64, folder: Option<i64>) -> Result<(), String> {
        self.conn.execute("UPDATE documents SET folder_id = ?2 WHERE id = ?1", params![id, folder]).map_err(|e| e.to_string())?;
        Ok(())
    }

    /// Deletes a folder; its documents and subfolders move to its parent.
    pub fn delete_folder(&mut self, id: i64) -> Result<(), String> {
        let tx = self.conn.transaction().map_err(|e| e.to_string())?;
        let parent: Option<i64> = tx.query_row("SELECT parent_id FROM folders WHERE id = ?1", params![id], |r| r.get(0)).map_err(|e| e.to_string())?;
        tx.execute("UPDATE documents SET folder_id = ?2 WHERE folder_id = ?1", params![id, parent]).map_err(|e| e.to_string())?;
        tx.execute("UPDATE folders SET parent_id = ?2 WHERE parent_id = ?1", params![id, parent]).map_err(|e| e.to_string())?;
        tx.execute("DELETE FROM folders WHERE id = ?1", params![id]).map_err(|e| e.to_string())?;
        tx.commit().map_err(|e| e.to_string())
    }

    pub fn setting(&self, key: &str) -> Option<String> {
        self.conn
            .query_row("SELECT value FROM settings WHERE key = ?1", params![key], |r| r.get(0))
            .optional()
            .ok()
            .flatten()
    }

    pub fn set_setting(&self, key: &str, value: &str) -> Result<(), String> {
        self.conn
            .execute(
                "INSERT INTO settings (key, value) VALUES (?1, ?2) ON CONFLICT(key) DO UPDATE SET value = excluded.value",
                params![key, value],
            )
            .map_err(|e| e.to_string())?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn store_and_load() {
        let path = std::env::temp_dir().join(format!("pcr_studio_db_test_{}.db", std::process::id()));
        let _ = std::fs::remove_file(&path);
        let mut db = Db::open(&path).unwrap();
        let mut r = Row::from_gapped("s1".into(), "desc".into(), b"--ACG-T".to_vec());
        r.set_meta("host", "human".into());
        let doc = Document::new(vec![r, Row::from_gapped("s2".into(), String::new(), b"AAACGTT".to_vec())]);
        let id = db.insert_document("test", DocKind::Alignment, &doc, Some(1), "", None, &|_| {}).unwrap();
        let (info, loaded, anns) = db.load_document(id, &|_| {}).unwrap();
        assert_eq!(info.reference, Some(1));
        assert_eq!(loaded.width, 7);
        assert_eq!(loaded.rows[0].start, 2);
        assert_eq!(loaded.rows[0].meta_value("host"), Some("human"));
        assert!(anns.is_empty());
        let aid = db
            .add_annotation(id, &Annotation { id: 0, name: "F1".into(), kind: AnnotationKind::ForwardPrimer, start: 1, end: 4, note: String::new() })
            .unwrap();
        assert_eq!(db.annotations(id).unwrap()[0].id, aid);
        let edited = Annotation { id: aid, name: "R1".into(), kind: AnnotationKind::ReversePrimer, start: 2, end: 6, note: "x".into() };
        db.update_annotation(&edited).unwrap();
        assert_eq!(db.annotations(id).unwrap(), vec![edited]);
        db.set_setting("k", "v1").unwrap();
        db.set_setting("k", "v2").unwrap();
        assert_eq!(db.setting("k").as_deref(), Some("v2"));
        // Folders: nesting, moving, cycle protection, deletion keeps contents.
        let a = db.create_folder("A", None).unwrap();
        let b = db.create_folder("B", Some(a)).unwrap();
        db.move_document(id, Some(b)).unwrap();
        assert_eq!(db.list_documents().unwrap()[0].folder, Some(b));
        assert!(db.move_folder(a, Some(b)).is_err());
        db.set_folder_collapsed(a, true).unwrap();
        assert!(db.list_folders().unwrap().iter().find(|f| f.id == a).unwrap().collapsed);
        db.delete_folder(b).unwrap();
        assert_eq!(db.list_documents().unwrap()[0].folder, Some(a));
        db.delete_folder(a).unwrap();
        assert_eq!(db.list_documents().unwrap()[0].folder, None);
        assert!(db.list_folders().unwrap().is_empty());
        // Editing is limited to single-sequence documents.
        assert!(db.replace_single_row(id, "x", DocKind::Sequences, &Row::default()).is_err());
        let one = Document::new(vec![Row::from_gapped("s".into(), String::new(), b"ACGTACGT".to_vec())]);
        let single = db.insert_document("single", DocKind::Sequences, &one, None, "", None, &|_| {}).unwrap();
        for (s, e) in [(0, 4), (2, 8)] {
            db.add_annotation(single, &Annotation { id: 0, name: "a".into(), kind: AnnotationKind::Region, start: s, end: e, note: String::new() })
                .unwrap();
        }
        let mut edited = Row::from_gapped("s2".into(), "d".into(), b"AC-GTA".to_vec());
        edited.set_meta("k", "v".into());
        assert_eq!(db.replace_single_row(single, "renamed", DocKind::Alignment, &edited).unwrap(), 1);
        let (info, loaded, anns) = db.load_document(single, &|_| {}).unwrap();
        assert_eq!((info.name.as_str(), info.kind, info.width, info.n_rows), ("renamed", DocKind::Alignment, 6, 1));
        assert_eq!((loaded.rows[0].name.as_str(), loaded.rows[0].data.as_slice()), ("s2", &b"AC-GTA"[..]));
        assert_eq!(loaded.rows[0].meta_value("k"), Some("v"));
        assert_eq!(anns.len(), 1);
        db.delete_document(single).unwrap();
        // Reordering rows: two cycles and a fixed row; metadata moves along.
        let names = ["a", "b", "c", "d", "e", "f"];
        let rows = names.iter().map(|n| {
            let mut r = Row::from_gapped(n.to_string(), String::new(), b"ACGT".to_vec());
            r.set_meta("n", n.to_string());
            r
        });
        let six = db.insert_document("six", DocKind::Sequences, &Document::new(rows.collect()), Some(2), "", None, &|_| {}).unwrap();
        db.reorder_rows(six, &[2, 0, 1, 3, 5, 4], Some(0)).unwrap();
        let (info, loaded, _) = db.load_document(six, &|_| {}).unwrap();
        let order: Vec<&str> = loaded.rows.iter().map(|r| r.name.as_str()).collect();
        assert_eq!(order, ["c", "a", "b", "d", "f", "e"]);
        assert!(loaded.rows.iter().all(|r| r.meta_value("n") == Some(r.name.as_str())));
        assert_eq!(info.reference, Some(0));
        assert!(db.reorder_rows(six, &[0, 1, 2], None).is_err());
        assert!(db.reorder_rows(six, &[0, 0, 1, 2, 3, 4], None).is_err());
        db.delete_document(six).unwrap();
        // Text documents: stored with their content; a missing folder becomes the top level.
        let text = db.insert_text_document("out", "line 1
line 2", "tool output", Some(9999)).unwrap();
        let info = db.document_info(text).unwrap().unwrap();
        assert_eq!((info.kind, info.n_rows, info.folder), (DocKind::Text, 0, None));
        assert_eq!(db.text_content(text).unwrap(), "line 1
line 2");
        assert_eq!(db.now_local().len(), 19);
        db.delete_document(text).unwrap();
        assert!(db.text_content(text).is_err());
        db.delete_document(id).unwrap();
        assert!(db.list_documents().unwrap().is_empty());
        assert!(db.annotations(id).unwrap().is_empty());
        drop(db);
        let _ = std::fs::remove_file(&path);
    }
}
