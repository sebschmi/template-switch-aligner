use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{
        alphabet::Alphabet,
        sequence::{GenomeSequence, OwnedGenomeSequence},
    },
};
use generic_a_star::{AStar, AStarResult, cost::AStarCost};
use indicatif::ProgressBar;
use lib_tsalign::a_star_aligner::{
    alignment_result::AlignmentResult,
    template_switch_distance::{EqualCostRange, TemplateSwitchDirection},
};
use log::{debug, trace};
use std::{
    fmt::Write,
    iter,
    time::{Duration, Instant},
};

use crate::{
    alignment::{
        Alignment, AlignmentType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
    },
    anchors::Anchors,
    chain_align::chainer::{Context, Identifier},
    chaining_cost_function::ChainingCostFunction,
    costs::AlignmentCosts,
    exact_chaining::{
        gap_affine::GapAffineAlignment, ts_12_jump::Ts12JumpAlignment,
        ts_34_jump::Ts34JumpAlignment,
    },
};

mod chainer;

#[expect(clippy::too_many_arguments)]
pub fn align<AlphabetType: Alphabet, Cost: AStarCost>(
    sequences: &AlignmentSequences,
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    alignment_costs: &AlignmentCosts<Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    max_match_run: u32,
    anchors: &Anchors,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) -> AlignmentResult<lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType, Cost> {
    let progress_bar = ProgressBar::new_spinner();
    progress_bar.enable_steady_tick(Duration::from_millis(200));

    let start_time = Instant::now();
    let mut chaining_duration = Duration::default();
    let mut evaluation_duration = Duration::default();

    let k = usize::try_from(max_match_run + 1).unwrap();
    let context = Context::new(anchors, chaining_cost_function);
    let mut astar = AStar::new(context);
    let mut chaining_execution_count = 0;
    let mut current_cost = Cost::zero();

    let (alignments, result) = loop {
        progress_bar.inc(1);
        progress_bar.set_message(format!(
            "Computing chain number {} (cost {})",
            chaining_execution_count + 1,
            current_cost,
        ));
        let chaining_start_time = Instant::now();

        astar.reset();
        astar.initialise();
        let result = astar.search();
        chaining_execution_count += 1;
        let chain = match result {
            AStarResult::FoundTarget { cost, .. } => {
                trace!("Found chain with cost {cost}");
                current_cost = cost;

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
                            Identifier::Secondary { index, ts_kind } => write!(
                                s,
                                "S{}{}",
                                ts_kind.digits(),
                                anchors.secondary(*ts_kind)[*index]
                            )
                            .unwrap(),
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

        let chaining_end_time = Instant::now();
        chaining_duration += chaining_end_time - chaining_start_time;

        let evaluation_start_time = Instant::now();

        let mut cost_increased = false;
        let mut alignments = Vec::new();
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
                    trace!("Aligning from start to end costs {}", alignment.cost());
                    cost_increased = astar
                        .context_mut()
                        .chaining_cost_function
                        .update_start_to_end(alignment.cost())
                        || cost_increased;
                    alignments.push(alignment.alignment().clone());
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
                    trace!(
                        "Aligning from start to P{index}{} costs {}",
                        anchors.primary[index],
                        alignment.cost()
                    );
                    cost_increased = astar
                        .context_mut()
                        .chaining_cost_function
                        .update_primary_from_start(index, alignment.cost())
                        || cost_increased;
                    alignments.push(alignment.alignment().clone());
                }
                (Identifier::Start, Identifier::Secondary { index, ts_kind }) => {
                    let end = anchors.secondary(ts_kind)[index].start(ts_kind);
                    let alignment = Ts12JumpAlignment::new(
                        start,
                        end,
                        sequences,
                        alignment_costs,
                        rc_fn,
                        max_match_run,
                    );
                    trace!(
                        "Aligning from start to S{}[{index}]{} costs {}",
                        ts_kind.digits(),
                        anchors.secondary(ts_kind)[index],
                        alignment.cost()
                    );
                    cost_increased = astar
                        .context_mut()
                        .chaining_cost_function
                        .update_jump_12_from_start(index, ts_kind, alignment.cost())
                        || cost_increased;
                    alignments.push(alignment.alignment().clone());
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
                    trace!(
                        "Aligning from P{index}{} to end costs {}",
                        anchors.primary[index],
                        alignment.cost()
                    );
                    cost_increased = astar
                        .context_mut()
                        .chaining_cost_function
                        .update_primary_to_end(index, alignment.cost())
                        || cost_increased;
                    alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                    alignments.push(alignment.alignment().clone());
                }
                (Identifier::Secondary { index, ts_kind }, Identifier::End) => {
                    let start = anchors.secondary(ts_kind)[index].end(ts_kind, k);
                    let alignment = Ts34JumpAlignment::new(
                        start,
                        end,
                        sequences,
                        alignment_costs,
                        rc_fn,
                        max_match_run,
                    );
                    trace!(
                        "Aligning from S{}[{index}]{} to end costs {}",
                        ts_kind.digits(),
                        anchors.secondary(ts_kind)[index],
                        alignment.cost()
                    );
                    cost_increased = astar
                        .context_mut()
                        .chaining_cost_function
                        .update_jump_34_to_end(index, ts_kind, alignment.cost())
                        || cost_increased;
                    alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                    alignments.push(alignment.alignment().clone());
                }
                (
                    Identifier::Primary { index: from_index },
                    Identifier::Primary { index: to_index },
                ) => {
                    if anchors.primary[from_index]
                        .is_direct_predecessor_of(&anchors.primary[to_index])
                    {
                        alignments.push(Alignment::from(vec![AlignmentType::Match]));
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
                    trace!(
                        "Aligning from P{from_index}{} to P{to_index}{} costs {}",
                        anchors.primary[from_index],
                        anchors.primary[to_index],
                        alignment.cost()
                    );
                    cost_increased = astar.context_mut().chaining_cost_function.update_primary(
                        from_index,
                        to_index,
                        alignment.cost(),
                    ) || cost_increased;
                    alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                    alignments.push(alignment.alignment().clone());
                }
                (
                    Identifier::Primary { index: from_index },
                    Identifier::Secondary {
                        index: to_index,
                        ts_kind,
                    },
                ) => {
                    let start = anchors.primary[from_index].end(k);
                    let end = anchors.secondary(ts_kind)[to_index].start(ts_kind);
                    let alignment = Ts12JumpAlignment::new(
                        start,
                        end,
                        sequences,
                        alignment_costs,
                        rc_fn,
                        max_match_run,
                    );
                    cost_increased = astar.context_mut().chaining_cost_function.update_jump_12(
                        from_index,
                        to_index,
                        ts_kind,
                        alignment.cost(),
                    ) || cost_increased;
                    trace!(
                        "Aligning from P{from_index}{} to S{}[{to_index}]{} costs {}",
                        anchors.primary[from_index],
                        ts_kind.digits(),
                        anchors.secondary(ts_kind)[to_index],
                        alignment.cost()
                    );
                    alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                    alignments.push(alignment.alignment().clone());
                }
                (
                    Identifier::Secondary {
                        index: from_index,
                        ts_kind,
                    },
                    Identifier::Secondary {
                        index: to_index,
                        ts_kind: to_ts_kind,
                    },
                ) => {
                    assert_eq!(ts_kind, to_ts_kind);
                    if anchors.secondary(ts_kind)[from_index]
                        .is_direct_predecessor_of(&anchors.secondary(ts_kind)[to_index])
                    {
                        alignments.push(Alignment::from(vec![AlignmentType::Match]));
                        continue 'update_chain;
                    }
                    let start = anchors.secondary(ts_kind)[from_index].end(ts_kind, k);
                    let end = anchors.secondary(ts_kind)[to_index].start(ts_kind);
                    let alignment = GapAffineAlignment::new(
                        start,
                        end,
                        sequences,
                        &alignment_costs.secondary_costs,
                        rc_fn,
                        max_match_run,
                    );
                    trace!(
                        "Aligning from S{}[{from_index}]{} to S{}[{to_index}]{} costs {}",
                        ts_kind.digits(),
                        anchors.secondary(ts_kind)[from_index],
                        ts_kind.digits(),
                        anchors.secondary(ts_kind)[to_index],
                        alignment.cost()
                    );
                    cost_increased = astar.context_mut().chaining_cost_function.update_secondary(
                        from_index,
                        to_index,
                        ts_kind,
                        alignment.cost(),
                    ) || cost_increased;
                    alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                    alignments.push(alignment.alignment().clone());
                }
                (
                    Identifier::Secondary {
                        index: from_index,
                        ts_kind,
                    },
                    Identifier::Primary { index: to_index },
                ) => {
                    let start = anchors.secondary(ts_kind)[from_index].end(ts_kind, k);
                    let end = anchors.primary[to_index].start();
                    let alignment = Ts34JumpAlignment::new(
                        start,
                        end,
                        sequences,
                        alignment_costs,
                        rc_fn,
                        max_match_run,
                    );
                    cost_increased = astar.context_mut().chaining_cost_function.update_jump_34(
                        from_index,
                        to_index,
                        ts_kind,
                        alignment.cost(),
                    ) || cost_increased;
                    trace!(
                        "Aligning from S{}[{from_index}]{} to P{to_index}{} (S{} to P{}) costs {}",
                        ts_kind.digits(),
                        anchors.secondary(ts_kind)[from_index],
                        anchors.primary[to_index],
                        start,
                        end,
                        alignment.cost()
                    );
                    alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                    alignments.push(alignment.alignment().clone());
                }
                (Identifier::End, _) | (_, Identifier::Start) => unreachable!(),
            }
        }

        let evaluation_end_time = Instant::now();
        evaluation_duration += evaluation_end_time - evaluation_start_time;

        if !cost_increased {
            break (alignments, result);
        }
    };

    progress_bar.finish_and_clear();
    debug!("Computed {chaining_execution_count} chains");
    debug!("Chaining took {:.1}s", chaining_duration.as_secs_f64());
    debug!("Evaluation took {:.1}s", evaluation_duration.as_secs_f64());

    let mut tsalign_alignment =
        lib_tsalign::a_star_aligner::alignment_result::alignment::Alignment::new();
    let mut is_primary = true;
    let mut anti_primary_gap = 0;

    for alignment in alignments {
        use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType as TsAlignAlignmentType;

        for (multiplicity, alignment_type) in alignment.alignment {
            match alignment_type {
                AlignmentType::Match => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimaryMatch
                    } else {
                        anti_primary_gap -= 1;
                        TsAlignAlignmentType::SecondaryMatch
                    },
                ),
                AlignmentType::Substitution => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimarySubstitution
                    } else {
                        anti_primary_gap -= 1;
                        TsAlignAlignmentType::SecondarySubstitution
                    },
                ),
                AlignmentType::GapA => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimaryInsertion
                    } else {
                        TsAlignAlignmentType::SecondaryInsertion
                    },
                ),
                AlignmentType::GapB => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimaryDeletion
                    } else {
                        anti_primary_gap -= 1;
                        TsAlignAlignmentType::SecondaryDeletion
                    },
                ),
                AlignmentType::TsStart { jump, ts_kind } => {
                    assert!(is_primary);
                    assert_eq!(multiplicity, 1);
                    is_primary = false;
                    anti_primary_gap = jump;

                    tsalign_alignment.push_n(
                        multiplicity,
                        TsAlignAlignmentType::TemplateSwitchEntrance {
                            first_offset: jump,
                            equal_cost_range: EqualCostRange::new_invalid(),
                            primary: ts_kind.descendant.into_tsalign_primary(),
                            secondary: ts_kind.ancestor.into_tsalign_secondary(),
                            direction: TemplateSwitchDirection::Reverse,
                        },
                    );
                }
                AlignmentType::TsEnd { jump } => {
                    assert!(!is_primary);
                    assert_eq!(multiplicity, 1);
                    is_primary = true;

                    tsalign_alignment.push_n(
                        multiplicity,
                        TsAlignAlignmentType::TemplateSwitchExit {
                            anti_primary_gap: anti_primary_gap + jump,
                        },
                    )
                }
            }
        }
    }

    let end_time = Instant::now();
    let duration_seconds = (end_time - start_time).as_secs_f64();

    AlignmentResult::new_with_target::<AlphabetType, _>(
        tsalign_alignment.into_inner(),
        VectorGenome::from_slice_u8(sequences.seq1())
            .unwrap()
            .as_genome_subsequence(),
        VectorGenome::from_slice_u8(sequences.seq2())
            .unwrap()
            .as_genome_subsequence(),
        sequences.seq1_name(),
        sequences.seq2_name(),
        start.primary_ordinate_a().unwrap(),
        start.primary_ordinate_b().unwrap(),
        result.without_node_identifier(),
        duration_seconds,
        0,
        0,
        0,
        sequences.seq1().len(),
        sequences.seq2().len(),
    )
}
