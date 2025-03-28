#[derive(Debug, Clone)]
pub struct LabelledSequence<'a> {
    label_base: String,
    labels: Vec<String>,
    sequence: &'a str,
    /// True if this sequence was created by taking a substring from another sequence with `limit < offset`.
    negative_substring: bool,
}
impl<'a> LabelledSequence<'a> {
    pub fn new(label_base: impl ToString, sequence: &'a str) -> Self {
        let label_base = label_base.to_string();
        let labels = vec![label_base.clone()];

        Self {
            label_base,
            labels,
            sequence,
            negative_substring: false,
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

    pub fn is_negative_substring(&self) -> bool {
        self.negative_substring
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

    /// Creates a substring from this string.
    ///
    /// If `offset > limit`, then the substring `limit..offset` is returned,
    /// and [`is_negative_substring()`](Self::is_negative_substring) will return `true` for the returned string.
    pub fn substring(&self, offset: usize, limit: usize) -> Self {
        assert!(!self.negative_substring);

        if offset <= limit {
            Self {
                label_base: self.label_base.clone(),
                labels: self.labels.clone(),
                sequence: char_substring(self.sequence, offset, limit),
                negative_substring: false,
            }
        } else {
            Self {
                label_base: self.label_base.clone(),
                labels: self.labels.clone(),
                sequence: char_substring(self.sequence, limit, offset),
                negative_substring: true,
            }
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
