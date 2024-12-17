use std::fmt::Display;

use super::{Chain, ChainLink};

impl Display for Chain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Chain [")?;
        for (index, link) in self.chain.iter().enumerate() {
            if index != 0 {
                write!(f, ", ")?;
            }

            write!(f, "{link}")?;
        }
        write!(f, "]")
    }
}

impl Display for ChainLink {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}; {})", self.identifier, self.cost)
    }
}
