use context::Context;
use generic_a_star::{AStar, cost::AStarCost};
use get_size2::GetSize;
use log::info;
use node::EdgeType;

use crate::seed::ChainingAnchors;

pub use context::ChainingCostsProvider;
pub use node::Identifier;

mod context;
mod display;
mod node;

#[derive(GetSize)]
pub struct Chain<Cost> {
    chain: Vec<ChainLink<Cost>>,
}

#[derive(GetSize)]
pub struct ChainLink<Cost> {
    identifier: Identifier,
    cost: Cost,
}

impl<Cost: AStarCost> Chain<Cost> {
    pub fn compute_chain<ChainingCosts: ChainingCostsProvider<Cost = Cost> + GetSize>(
        chaining_costs: ChainingCosts,
        chaining_anchors: ChainingAnchors,
    ) -> Self {
        info!("Computing chain...");
        let mut a_star = AStar::new(Context::new(chaining_costs, chaining_anchors));
        a_star.initialise();
        a_star.search();

        let mut chain = Vec::new();
        chain.extend(
            a_star
                .backtrack_with_costs()
                .map(|(EdgeType { to, .. }, cost)| match to {
                    Identifier::Root => unreachable!(),
                    identifier => ChainLink { identifier, cost },
                }),
        );
        chain.push(ChainLink {
            identifier: Identifier::Root,
            cost: Cost::zero(),
        });
        chain.reverse();

        // Transform costs from cost to root into cost to target.
        let total_cost = chain.last().unwrap().cost;
        chain
            .iter_mut()
            .for_each(|link| link.cost = total_cost - link.cost);

        // Assert that costs and blocks are ordered, and blocks do not overlap.
        debug_assert!({
            chain.windows(2).all(|window| {
                let first = &window[0];
                let second = &window[1];

                match (&first.identifier, &second.identifier) {
                    (Identifier::Root, Identifier::Root)
                    | (Identifier::Anchor { .. }, Identifier::Root)
                    | (Identifier::Target, Identifier::Root)
                    | (Identifier::Target, Identifier::Anchor { .. })
                    | (Identifier::Target, Identifier::Target) => false,
                    (Identifier::Root, Identifier::Anchor { .. })
                    | (Identifier::Root, Identifier::Target)
                    | (Identifier::Anchor { .. }, Identifier::Target) => first.cost >= second.cost,
                    (
                        Identifier::Anchor {
                            anchor: first_anchor,
                        },
                        Identifier::Anchor {
                            anchor: second_anchor,
                        },
                    ) => {
                        first.cost <= second.cost
                            && first_anchor.reference_block().end
                                <= second_anchor.reference_block().start
                            && first_anchor.query_block().end <= second_anchor.query_block().start
                    }
                }
            })
        });

        Self { chain }
    }

    pub fn chain_lower_bound(&self, reference_index: usize, query_index: usize) -> Cost {
        match self.chain.binary_search_by_key(
            &(reference_index, query_index),
            |ChainLink { identifier, .. }| match identifier {
                Identifier::Root => (0, 0),
                Identifier::Anchor { anchor } => (
                    anchor.reference_block().end - 1,
                    anchor.query_block().end - 1,
                ),
                Identifier::Target => (usize::MAX, usize::MAX),
            },
        ) {
            Ok(index) => self.chain[index].cost,
            Err(index) => self
                .chain
                .get(index)
                .map(|link| link.cost)
                .unwrap_or(Cost::zero()),
        }
    }
}
