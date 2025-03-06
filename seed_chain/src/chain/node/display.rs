use std::fmt::Display;

use super::{Identifier, Node};

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}; ", self.identifier)?;
        if let Some(predecessor) = &self.predecessor {
            write!(f, "{predecessor}")?;
        } else {
            write!(f, "root")?;
        }
        write!(f, "; {}]", self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Identifier::Root => write!(f, "[root]"),
            Identifier::Anchor { anchor } => write!(
                f,
                "[R: {}..{}; Q: {}..{}]",
                anchor.reference_block().start,
                anchor.reference_block().end,
                anchor.query_block().start,
                anchor.query_block().end
            ),
            Identifier::Target => write!(f, "[target]"),
        }
    }
}
