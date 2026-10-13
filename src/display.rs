//! Display list: which rows are shown, in which order, grouped and collapsed.

use crate::model::Document;
use crate::seq::GAP;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SortKey {
    Original,
    Name,
    StartPosition,
    Length,
    IdentityToReference,
    Meta(String),
}

impl std::fmt::Display for SortKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SortKey::Original => f.write_str("Document order"),
            SortKey::Name => f.write_str("Name"),
            SortKey::StartPosition => f.write_str("Start position"),
            SortKey::Length => f.write_str("Sequence length"),
            SortKey::IdentityToReference => f.write_str("Identity to reference"),
            SortKey::Meta(k) => write!(f, "Metadata: {k}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DisplayOptions {
    pub sort: SortKey,
    pub descending: bool,
    pub group_by: Option<String>,
    pub collapse_identical: bool,
    pub name_filter: String,
    pub folded_groups: HashSet<String>,
}

impl Default for DisplayOptions {
    fn default() -> Self {
        DisplayOptions {
            sort: SortKey::Original,
            descending: false,
            group_by: None,
            collapse_identical: false,
            name_filter: String::new(),
            folded_groups: HashSet::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Item {
    Group { label: String, sequences: usize, folded: bool },
    /// A sequence row; `count` > 1 when identical sequences are collapsed into it.
    Seq { row: usize, count: usize },
}

/// Compares strings numerically when both parse as numbers, otherwise
/// case-insensitively with runs of digits compared by value ("seq9" < "seq10").
pub fn natural_cmp(a: &str, b: &str) -> Ordering {
    match (a.trim().parse::<f64>(), b.trim().parse::<f64>()) {
        (Ok(x), Ok(y)) => x.partial_cmp(&y).unwrap_or(Ordering::Equal),
        (Ok(_), Err(_)) => Ordering::Less,
        (Err(_), Ok(_)) => Ordering::Greater,
        _ => digit_runs_cmp(a.to_lowercase().as_bytes(), b.to_lowercase().as_bytes()),
    }
}

fn digit_runs_cmp(mut x: &[u8], mut y: &[u8]) -> Ordering {
    loop {
        match (x.first(), y.first()) {
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(p), Some(q)) if p.is_ascii_digit() && q.is_ascii_digit() => {
                let n = x.iter().take_while(|c| c.is_ascii_digit()).count();
                let m = y.iter().take_while(|c| c.is_ascii_digit()).count();
                let trim = |d: &[u8]| -> usize { d.iter().take_while(|&&c| c == b'0').count() };
                let (dx, dy) = (&x[trim(&x[..n])..n], &y[trim(&y[..m])..m]);
                let o = dx.len().cmp(&dy.len()).then_with(|| dx.cmp(dy));
                if o != Ordering::Equal {
                    return o;
                }
                (x, y) = (&x[n..], &y[m..]);
            }
            (Some(p), Some(q)) => {
                if p != q {
                    return p.cmp(q);
                }
                (x, y) = (&x[1..], &y[1..]);
            }
        }
    }
}

/// Permanent sort of a document's sequences (right-click menu).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RowSort {
    Name,
    /// Ungapped length; longest first unless reversed.
    Length,
    /// Bases in these alignment columns, compared left to right.
    Bases(Vec<usize>),
}

impl RowSort {
    /// Direction used the first time: names A→Z, longest first, bases A→T.
    pub fn default_descending(&self) -> bool {
        matches!(self, RowSort::Length)
    }
}

/// Sort rank of an alignment cell: A, C, G, T, other IUPAC codes
/// (alphabetically), gap, no coverage.
pub fn base_rank(b: Option<u8>) -> u8 {
    match b.map(|b| b.to_ascii_uppercase()) {
        Some(b'A') => 0,
        Some(b'C') => 1,
        Some(b'G') => 2,
        Some(b'T' | b'U') => 3,
        Some(GAP) => 101,
        Some(b) if b.is_ascii_uppercase() => 4 + (b - b'A'),
        Some(_) => 100,
        None => 102,
    }
}

/// New row order for sorting a document by `key`: row `order[i]` becomes
/// row `i`. All sequences are sorted; the reference goes first. The sort is
/// stable, so ties keep their current order.
pub fn sort_order(doc: &Document, reference: Option<usize>, key: &RowSort, descending: bool) -> Vec<usize> {
    let mut order: Vec<usize> = (0..doc.rows.len()).filter(|&i| Some(i) != reference).collect();
    let lens: Vec<usize> = match key {
        RowSort::Length => doc.rows.iter().map(|r| r.ungapped_len()).collect(),
        _ => Vec::new(),
    };
    order.sort_by(|&a, &b| {
        let (ra, rb) = (&doc.rows[a], &doc.rows[b]);
        let o = match key {
            RowSort::Name => natural_cmp(&ra.name, &rb.name),
            RowSort::Length => lens[a].cmp(&lens[b]),
            RowSort::Bases(cols) => cols.iter().map(|&c| base_rank(ra.at(c))).cmp(cols.iter().map(|&c| base_rank(rb.at(c)))),
        };
        if descending { o.reverse() } else { o }
    });
    if let Some(r) = reference {
        order.insert(0, r);
    }
    order
}

/// Builds the display list (the reference row is shown pinned and excluded here).
/// `identity` holds per-row identity to the reference, used for sorting.
pub fn build(doc: &Document, reference: Option<usize>, opts: &DisplayOptions, identity: Option<&[f32]>) -> Vec<Item> {
    let filter = opts.name_filter.trim().to_lowercase();
    let mut rows: Vec<usize> = (0..doc.rows.len())
        .filter(|&i| Some(i) != reference)
        .filter(|&i| filter.is_empty() || doc.rows[i].name.to_lowercase().contains(&filter) || doc.rows[i].description.to_lowercase().contains(&filter))
        .collect();

    let cmp = |a: &usize, b: &usize| -> Ordering {
        let (ra, rb) = (&doc.rows[*a], &doc.rows[*b]);
        let o = match &opts.sort {
            SortKey::Original => a.cmp(b),
            SortKey::Name => natural_cmp(&ra.name, &rb.name),
            SortKey::StartPosition => ra.start.cmp(&rb.start).then(ra.end().cmp(&rb.end())),
            SortKey::Length => ra.ungapped_len().cmp(&rb.ungapped_len()),
            SortKey::IdentityToReference => match identity {
                Some(id) => id[*a].partial_cmp(&id[*b]).unwrap_or(Ordering::Equal),
                None => Ordering::Equal,
            },
            SortKey::Meta(k) => natural_cmp(ra.meta_value(k).unwrap_or(""), rb.meta_value(k).unwrap_or("")),
        };
        let o = if opts.descending { o.reverse() } else { o };
        o.then(a.cmp(b))
    };
    rows.sort_by(cmp);

    let collapse = |members: &[usize]| -> Vec<Item> {
        if !opts.collapse_identical {
            return members.iter().map(|&row| Item::Seq { row, count: 1 }).collect();
        }
        let mut seen: HashMap<(usize, &[u8]), usize> = HashMap::new();
        let mut out: Vec<Item> = Vec::new();
        for &row in members {
            let r = &doc.rows[row];
            match seen.get(&(r.start, r.data.as_slice())) {
                Some(&pos) => {
                    if let Item::Seq { count, .. } = &mut out[pos] {
                        *count += 1;
                    }
                }
                None => {
                    seen.insert((r.start, r.data.as_slice()), out.len());
                    out.push(Item::Seq { row, count: 1 });
                }
            }
        }
        out
    };

    match &opts.group_by {
        None => collapse(&rows),
        Some(key) => {
            let mut order: Vec<String> = Vec::new();
            let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
            for &r in &rows {
                let label = doc.rows[r].meta_value(key).unwrap_or("(none)").to_string();
                groups.entry(label.clone()).or_insert_with(|| {
                    order.push(label.clone());
                    Vec::new()
                });
                groups.get_mut(&label).unwrap().push(r);
            }
            order.sort_by(|a, b| natural_cmp(a, b));
            let mut out = Vec::new();
            for label in order {
                let members = &groups[&label];
                let folded = opts.folded_groups.contains(&label);
                out.push(Item::Group { label: label.clone(), sequences: members.len(), folded });
                if !folded {
                    out.extend(collapse(members));
                }
            }
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Row;

    fn doc() -> Document {
        let mk = |name: &str, s: &str, g: &str| {
            let mut r = Row::from_gapped(name.into(), String::new(), s.as_bytes().to_vec());
            r.set_meta("g", g.into());
            r
        };
        Document::new(vec![mk("ref", "ACGT", "-"), mk("b10", "ACGT", "x"), mk("b9", "ACGA", "y"), mk("a", "ACGT", "x")])
    }

    #[test]
    fn sort_collapse_group() {
        let d = doc();
        let mut o = DisplayOptions { sort: SortKey::Name, ..Default::default() };
        let items = build(&d, Some(0), &o, None);
        // Natural order: "b9" before "b10".
        assert_eq!(items, vec![Item::Seq { row: 3, count: 1 }, Item::Seq { row: 2, count: 1 }, Item::Seq { row: 1, count: 1 }]);
        o.collapse_identical = true;
        let items = build(&d, Some(0), &o, None);
        assert_eq!(items.len(), 2);
        assert_eq!(items[0], Item::Seq { row: 3, count: 2 });
        o.group_by = Some("g".into());
        o.collapse_identical = false;
        let items = build(&d, Some(0), &o, None);
        assert_eq!(items[0], Item::Group { label: "x".into(), sequences: 2, folded: false });
        assert_eq!(items.len(), 5);
        o.name_filter = "b".into();
        let items = build(&d, Some(0), &o, None);
        assert_eq!(items.len(), 4);
    }

    #[test]
    fn natural_order() {
        assert_eq!(natural_cmp("9", "10"), Ordering::Less);
        assert_eq!(natural_cmp("1.5", "1.25"), Ordering::Greater);
        assert_eq!(natural_cmp("b", "A"), Ordering::Greater);
        assert_eq!(natural_cmp("seq_9", "seq_10"), Ordering::Less);
        assert_eq!(natural_cmp("Seq_010", "seq_10"), Ordering::Equal);
        assert_eq!(natural_cmp("seq_10a", "seq_10"), Ordering::Greater);
        assert_eq!(natural_cmp("seq_10", "seq_1x"), Ordering::Greater);
    }

    #[test]
    fn sorting_rows() {
        let mk = |name: &str, s: &str| Row::from_gapped(name.into(), String::new(), s.as_bytes().to_vec());
        // Row 0 is the reference; leading/trailing gaps are no coverage.
        let d = Document::new(vec![
            mk("ref", "ACGTAC"),
            mk("s10", "ACGTA-"),
            mk("s9", "--GTAC"),
            mk("S2", "ACNTAC"),
            mk("s1", "AC-TAC"),
            mk("s3", "ACTTAC"),
        ]);
        let names = |o: &[usize]| o.iter().map(|&i| d.rows[i].name.as_str()).collect::<Vec<_>>();
        assert_eq!(names(&sort_order(&d, Some(0), &RowSort::Name, false)), ["ref", "s1", "S2", "s3", "s9", "s10"]);
        assert_eq!(names(&sort_order(&d, Some(0), &RowSort::Name, true)), ["ref", "s10", "s9", "s3", "S2", "s1"]);
        // Without a reference every row is sorted.
        assert_eq!(names(&sort_order(&d, None, &RowSort::Name, false)), ["ref", "s1", "S2", "s3", "s9", "s10"]);
        // Longest first; ties keep the current order.
        assert_eq!(names(&sort_order(&d, Some(0), &RowSort::Length, true)), ["ref", "S2", "s3", "s10", "s1", "s9"]);
        assert_eq!(names(&sort_order(&d, Some(0), &RowSort::Length, false)), ["ref", "s9", "s10", "s1", "S2", "s3"]);
        // Column 2: G, G, N, gap, T → A C G T, others, gap, no coverage.
        let by2 = RowSort::Bases(vec![2]);
        assert_eq!(names(&sort_order(&d, Some(0), &by2, false)), ["ref", "s10", "s9", "s3", "S2", "s1"]);
        assert_eq!(names(&sort_order(&d, Some(0), &by2, true)), ["ref", "s1", "S2", "s3", "s10", "s9"]);
        // Column 0: s9 has no coverage and goes last.
        assert_eq!(names(&sort_order(&d, Some(0), &RowSort::Bases(vec![0]), false)), ["ref", "s10", "S2", "s1", "s3", "s9"]);
        // Several columns compare left to right.
        let sel = RowSort::Bases(vec![1, 2, 5]);
        assert_eq!(names(&sort_order(&d, Some(0), &sel, false)), ["ref", "s10", "s3", "S2", "s1", "s9"]);
        assert!(base_rank(Some(b'T')) < base_rank(Some(b'N')));
        assert!(base_rank(Some(b'N')) < base_rank(Some(GAP)));
        assert!(base_rank(Some(GAP)) < base_rank(None));
        assert_eq!(base_rank(Some(b'u')), base_rank(Some(b'T')));
    }
}
