use std::fmt::Display;

use super::{Chain, ChainLink};

impl<Cost: Display> Display for Chain<Cost> {
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

impl<Cost: Display> Display for ChainLink<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}; {})", self.identifier, self.cost)
    }
}
