#[derive(Debug, Clone)]
pub struct LabelledSequence<'a> {
    label_base: String,
    labels: Vec<String>,
    sequence: &'a str,
}
impl<'a> LabelledSequence<'a> {
    pub fn new(label_base: impl ToString, sequence: &'a str) -> Self {
        let label_base = label_base.to_string();
        let labels = vec![label_base.clone()];

        Self {
            label_base,
            labels,
            sequence,
        }
    }

    pub fn label(&self) -> &String {
        self.labels.last().unwrap()
    }

    pub fn label_cloned(&self) -> String {
        self.label().to_owned()
    }

    pub fn all_labels(&self) -> &[String] {
        &self.labels
    }

    pub fn sequence(&self) -> &'a str {
        self.sequence
    }

    /// Adds a new label created with the same base.
    ///
    /// If the current label was renamed by this operation, returns the renaming as `Some(old_name, new_name)`.
    pub fn increment_label(&mut self) -> Option<(String, String)> {
        let result = if self.labels.len() == 1 {
            let old_name = self.labels[0].clone();
            self.labels[0].push_str(" 1");
            let new_name = self.labels[0].clone();
            Some((old_name, new_name))
        } else {
            None
        };

        self.labels
            .push(format!("{} {}", self.label_base, self.labels.len() + 1));
        result
    }

    pub fn substring(&self, offset: usize, limit: usize) -> Self {
        Self {
            label_base: self.label_base.clone(),
            labels: self.labels.clone(),
            sequence: char_substring(self.sequence, offset, limit),
        }
    }
}

fn char_substring(string: &str, offset: usize, limit: usize) -> &str {
    //trace!("Taking substring {offset}..{limit} of {string}");

    let mut indices = string
        .char_indices()
        .map(|(index, _)| index)
        .chain(Some(string.len()));
    let byte_offset = indices.by_ref().nth(offset).unwrap_or_else(|| {
        panic!(
            "The string contains {} characters, but the offset is {offset}",
            string.chars().count()
        )
    });

    if offset == limit {
        &string[byte_offset..byte_offset]
    } else {
        assert!(offset < limit, "offset < limit: {offset} < {limit}");
        let byte_limit = indices.nth(limit - offset - 1).unwrap_or_else(|| {
            panic!(
                "The string contains {} characters, but the limit is {limit}",
                string.chars().count()
            )
        });
        &string[byte_offset..byte_limit]
    }
}
