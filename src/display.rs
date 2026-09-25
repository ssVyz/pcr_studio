//! Display list: which rows are shown, in which order, grouped and collapsed.

use crate::model::Document;
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
            SortKey::Original => f.write_str("Original order"),
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

/// Compares strings numerically when both parse as numbers.
pub fn natural_cmp(a: &str, b: &str) -> Ordering {
    match (a.trim().parse::<f64>(), b.trim().parse::<f64>()) {
        (Ok(x), Ok(y)) => x.partial_cmp(&y).unwrap_or(Ordering::Equal),
        (Ok(_), Err(_)) => Ordering::Less,
        (Err(_), Ok(_)) => Ordering::Greater,
        _ => a.to_lowercase().cmp(&b.to_lowercase()),
    }
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
        assert_eq!(items, vec![Item::Seq { row: 3, count: 1 }, Item::Seq { row: 1, count: 1 }, Item::Seq { row: 2, count: 1 }]);
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
        assert_eq!(natural_cmp("b", "A"), Ordering::Greater);
    }
}
