//! SQLite document library, stored next to the executable (`pcr_studio.db`).

use crate::model::{Annotation, AnnotationKind, DocInfo, DocKind, Document, Row};
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
CREATE TABLE IF NOT EXISTS settings (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
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
        Ok(Db { conn })
    }

    pub fn list_documents(&self) -> Result<Vec<DocInfo>, String> {
        let mut st = self
            .conn
            .prepare("SELECT id, name, kind, n_rows, width, reference, created, description FROM documents ORDER BY id")
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
        progress: &dyn Fn(f32),
    ) -> Result<i64, String> {
        let tx = self.conn.transaction().map_err(|e| e.to_string())?;
        tx.execute(
            "INSERT INTO documents (name, kind, n_rows, width, reference, description) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![name, kind.as_str(), doc.rows.len() as i64, doc.width as i64, reference.map(|r| r as i64), description],
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

    pub fn delete_annotation(&self, id: i64) -> Result<(), String> {
        self.conn.execute("DELETE FROM annotations WHERE id = ?1", params![id]).map_err(|e| e.to_string())?;
        Ok(())
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
        let id = db.insert_document("test", DocKind::Alignment, &doc, Some(1), "", &|_| {}).unwrap();
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
        db.set_setting("k", "v1").unwrap();
        db.set_setting("k", "v2").unwrap();
        assert_eq!(db.setting("k").as_deref(), Some("v2"));
        db.delete_document(id).unwrap();
        assert!(db.list_documents().unwrap().is_empty());
        assert!(db.annotations(id).unwrap().is_empty());
        drop(db);
        let _ = std::fs::remove_file(&path);
    }
}
