use generic_a_star::{AStar, AStarResult, cost::AStarCost};
use log::debug;
use std::fmt::Write;

use crate::{
    alignment::{coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    anchors::Anchors,
    chain_align::chainer::{Context, Identifier},
    chaining_cost_function::ChainingCostFunction,
};

mod chainer;

pub fn align<Cost: AStarCost>(
    sequences: &AlignmentSequences,
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    anchors: &Anchors,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) {
    let context = Context::new(anchors, chaining_cost_function);
    let mut astar = AStar::new(context);
    astar.initialise();

    loop {
        let chain = match astar.search() {
            AStarResult::FoundTarget { cost, .. } => {
                debug!("Found chain with cost {cost}");
                let mut chain = astar.reconstruct_path();
                chain.push(Identifier::End);
                astar.reset();
                debug!("Chain (len: {}):\n{}", chain.len(), {
                    let mut s = String::new();
                    let mut once = true;
                    for identifier in &chain {
                        if once {
                            once = false;
                        } else {
                            writeln!(s).unwrap();
                        }
                        match identifier {
                            Identifier::Start => write!(s, "start").unwrap(),
                            Identifier::Primary { index } => {
                                write!(s, "P{}", anchors.primary[*index]).unwrap()
                            }
                            Identifier::End => write!(s, "end").unwrap(),
                        }
                    }
                    s
                });
                chain
            }
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Memory limit is None"),
            AStarResult::NoTarget => panic!("No chain found"),
        };

        for window in chain.windows(2) {
            let from_anchor = window[0];
            let to_anchor = window[1];

            match (from_anchor, to_anchor) {
                (Identifier::Start, Identifier::End) => todo!(),
                (Identifier::Start, Identifier::Primary { index }) => todo!(),
                (Identifier::Primary { index }, Identifier::End) => todo!(),
                (
                    Identifier::Primary { index: from_index },
                    Identifier::Primary { index: to_index },
                ) => todo!(),
                (Identifier::End, _) | (_, Identifier::Start) => unreachable!(),
            }
        }

        todo!()
    }
}
