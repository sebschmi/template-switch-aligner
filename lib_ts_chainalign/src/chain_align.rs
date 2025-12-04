use generic_a_star::{AStar, AStarResult, cost::AStarCost};
use log::{debug, trace};
use std::fmt::Write;

use crate::{
    alignment::{coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    anchors::Anchors,
    chain_align::chainer::{Context, Identifier},
    chaining_cost_function::ChainingCostFunction,
    costs::AlignmentCosts,
    exact_chaining::gap_affine::GapAffineAlignment,
};

mod chainer;

#[expect(clippy::too_many_arguments)]
pub fn align<Cost: AStarCost>(
    sequences: &AlignmentSequences,
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    alignment_costs: &AlignmentCosts<Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    max_match_run: u32,
    anchors: &Anchors,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) {
    let k = usize::try_from(max_match_run + 1).unwrap();
    let context = Context::new(anchors, chaining_cost_function);
    let mut astar = AStar::new(context);

    loop {
        astar.reset();
        astar.initialise();
        let chain = match astar.search() {
            AStarResult::FoundTarget { cost, .. } => {
                debug!("Found chain with cost {cost}");
                let mut chain = astar.reconstruct_path();
                chain.push(Identifier::End);
                trace!("Chain (len: {}):\n{}", chain.len(), {
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

        let mut cost_increased = false;
        'update_chain: for window in chain.windows(2) {
            let from_anchor = window[0];
            let to_anchor = window[1];

            match (from_anchor, to_anchor) {
                (Identifier::Start, Identifier::End) => {
                    let alignment = GapAffineAlignment::new(
                        start,
                        end,
                        sequences,
                        &alignment_costs.primary_costs,
                        rc_fn,
                        max_match_run,
                    );
                    cost_increased = cost_increased
                        || astar
                            .context_mut()
                            .chaining_cost_function
                            .update_start_to_end(alignment.cost());
                }
                (Identifier::Start, Identifier::Primary { index }) => {
                    let end = anchors.primary[index].start();
                    let alignment = GapAffineAlignment::new(
                        start,
                        end,
                        sequences,
                        &alignment_costs.primary_costs,
                        rc_fn,
                        max_match_run,
                    );
                    cost_increased = cost_increased
                        || astar
                            .context_mut()
                            .chaining_cost_function
                            .update_primary_from_start(index, alignment.cost());
                }
                (Identifier::Primary { index }, Identifier::End) => {
                    let start = anchors.primary[index].end(k);
                    let alignment = GapAffineAlignment::new(
                        start,
                        end,
                        sequences,
                        &alignment_costs.primary_costs,
                        rc_fn,
                        max_match_run,
                    );
                    cost_increased = cost_increased
                        || astar
                            .context_mut()
                            .chaining_cost_function
                            .update_primary_to_end(index, alignment.cost());
                }
                (
                    Identifier::Primary { index: from_index },
                    Identifier::Primary { index: to_index },
                ) => {
                    if anchors.primary[from_index]
                        .is_direct_predecessor_of(&anchors.primary[to_index])
                    {
                        continue 'update_chain;
                    }

                    let start = anchors.primary[from_index].end(k);
                    let end = anchors.primary[to_index].start();
                    let alignment = GapAffineAlignment::new(
                        start,
                        end,
                        sequences,
                        &alignment_costs.primary_costs,
                        rc_fn,
                        max_match_run,
                    );
                    cost_increased = cost_increased
                        || astar.context_mut().chaining_cost_function.update_primary(
                            from_index,
                            to_index,
                            alignment.cost(),
                        );
                }
                (Identifier::End, _) | (_, Identifier::Start) => unreachable!(),
            }
        }

        if !cost_increased {
            break;
        }
    }

    todo!("alignment found")
}
