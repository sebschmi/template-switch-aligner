use std::hash::{BuildHasher, DefaultHasher};

/// A [`BuildHasher`] that builds the standard libraries default hasher with a deterministic seed.
///
/// See also [`DefaultHasher::new`].
#[derive(Default)]
pub struct DeterministicDefaultHasher;

impl BuildHasher for DeterministicDefaultHasher {
    type Hasher = DefaultHasher;

    fn build_hasher(&self) -> Self::Hasher {
        DefaultHasher::new()
    }
}
