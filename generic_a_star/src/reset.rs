use rustc_hash::FxHashMapSeed;

pub trait Reset {
    fn reset(&mut self);
}

impl Reset for () {
    fn reset(&mut self) {}
}

impl<Key, Value> Reset for FxHashMapSeed<Key, Value> {
    fn reset(&mut self) {
        self.clear();
    }
}
