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

impl<const N: usize, T: Reset> Reset for [T; N] {
    fn reset(&mut self) {
        self.iter_mut().for_each(Reset::reset);
    }
}
