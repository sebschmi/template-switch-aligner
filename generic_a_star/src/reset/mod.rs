use std::collections::HashMap;

pub trait Reset {
    fn reset(&mut self);
}

impl Reset for () {
    fn reset(&mut self) {}
}

impl<Key, Value, Hasher> Reset for HashMap<Key, Value, Hasher> {
    fn reset(&mut self) {
        self.clear();
    }
}
