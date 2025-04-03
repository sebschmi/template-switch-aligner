pub struct IndexedStr<'a> {
    index: Vec<usize>,
    string: &'a str,
}

impl<'a> IndexedStr<'a> {
    pub fn new(string: &'a str) -> Self {
        let index = string.char_indices().map(|(index, _)| index).collect();

        Self { index, string }
    }

    pub fn char_at(&self, index: usize) -> char {
        self.string[self.index[index]..].chars().next().unwrap()
    }
}
