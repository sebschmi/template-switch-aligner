pub trait Reset {
    fn reset(&mut self);
}

impl Reset for () {
    fn reset(&mut self) {}
}

impl<T> Reset for Option<T> {
    fn reset(&mut self) {
        *self = None;
    }
}

impl<A: Reset, B: Reset> Reset for (A, B) {
    fn reset(&mut self) {
        self.0.reset();
        self.1.reset();
    }
}
