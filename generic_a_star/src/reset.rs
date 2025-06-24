use std::collections::HashMap;

use memtally::Tracked;

pub trait Reset {
    fn reset(&mut self);
}

impl Reset for () {
    fn reset(&mut self) {}
}

impl<Key, Value, Hasher> Reset for Tracked<HashMap<Key, Value, Hasher>> {
    fn reset(&mut self) {
        self.clear();
    }
}
