use std::hash::{BuildHasher, DefaultHasher};

#[derive(Default)]
pub struct DeterministicDefaultHasher;

impl BuildHasher for DeterministicDefaultHasher {
    type Hasher = DefaultHasher;

    fn build_hasher(&self) -> Self::Hasher {
        DefaultHasher::new()
    }
}
