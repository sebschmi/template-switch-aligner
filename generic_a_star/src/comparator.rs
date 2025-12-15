use std::cmp::Ordering;

use compare::Compare;

use crate::AStarNode;

#[derive(Debug, Default)]
pub struct AStarNodeComparator;

impl<T: AStarNode> Compare<T> for AStarNodeComparator {
    fn compare(&self, l: &T, r: &T) -> Ordering {
        l.cmp(r).reverse().then_with(|| {
            l.secondary_maximisable_score()
                .cmp(&r.secondary_maximisable_score())
        })
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Display;

    use compare::Compare;

    use crate::{AStarIdentifier, AStarNode, cost::U64Cost};

    use super::AStarNodeComparator;

    #[derive(Debug, PartialEq, Eq)]
    struct Node {
        cost: U64Cost,
        lower_bound: U64Cost,
        secondary_maximisable_score: usize,
    }

    impl Node {
        fn new(cost: u64, lower_bound: u64, secondary_maximisable_score: usize) -> Self {
            Self {
                cost: cost.into(),
                lower_bound: lower_bound.into(),
                secondary_maximisable_score,
            }
        }
    }

    impl AStarNode for Node {
        type Identifier = ();

        type EdgeType = ();

        type Cost = U64Cost;

        fn identifier(&self) -> &Self::Identifier {
            unimplemented!()
        }

        fn cost(&self) -> Self::Cost {
            self.cost
        }

        fn a_star_lower_bound(&self) -> Self::Cost {
            self.lower_bound
        }

        fn secondary_maximisable_score(&self) -> usize {
            self.secondary_maximisable_score
        }

        fn predecessor(&self) -> Option<&Self::Identifier> {
            unimplemented!()
        }

        fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
            unimplemented!()
        }

        fn required_memory() -> usize {
            unimplemented!()
        }
    }

    impl Display for Node {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "{} + {}; {}",
                self.cost, self.lower_bound, self.secondary_maximisable_score
            )
        }
    }

    impl Ord for Node {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            (self.cost + self.lower_bound).cmp(&(other.cost + other.lower_bound))
        }
    }

    impl PartialOrd for Node {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    impl AStarIdentifier for () {}

    #[test]
    fn compare() {
        // Heap is a max heap, hence smaller nodes need to be bigger.
        assert!(AStarNodeComparator.compares_eq(&Node::new(4, 5, 6), &Node::new(4, 5, 6)));
        assert!(AStarNodeComparator.compares_gt(&Node::new(4, 4, 6), &Node::new(4, 5, 6)));
        assert!(AStarNodeComparator.compares_gt(&Node::new(4, 4, 6), &Node::new(4, 5, 6)));
        assert!(AStarNodeComparator.compares_gt(&Node::new(4, 5, 7), &Node::new(4, 5, 6)));
        assert!(AStarNodeComparator.compares_lt(&Node::new(4, 5, 5), &Node::new(4, 5, 6)));
        assert!(AStarNodeComparator.compares_lt(&Node::new(4, 6, 6), &Node::new(4, 5, 6)));
        assert!(AStarNodeComparator.compares_lt(&Node::new(5, 5, 6), &Node::new(4, 5, 6)));
    }
}
