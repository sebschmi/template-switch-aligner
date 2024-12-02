use std::{
    collections::{HashMap, HashSet},
    fmt::{Display, Write},
    iter,
};

use compact_genome::{
    implementation::vec_sequence::{SliceSubGenome, VectorGenome},
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use generic_a_star::{cost::Cost, AStar, AStarNode, AStarResult};
use log::{debug, trace};
use ndarray::Array2;

use crate::{
    a_star_aligner::{
        alignment_result::IAlignmentType,
        template_switch_distance::{
            strategies::{
                chaining::NoChainingStrategy, node_ord::CostOnlyNodeOrdStrategy,
                secondary_deletion_strategy::ForbidSecondaryDeletionStrategy,
                shortcut::NoShortcutStrategy,
                template_switch_count::MaxTemplateSwitchCountStrategy,
                template_switch_min_length::NoTemplateSwitchMinLengthStrategy,
                AlignmentStrategySelection,
            },
            Context, Identifier, Node,
        },
        AlignmentContext,
    },
    config::TemplateSwitchConfig,
    costs::gap_affine::GapAffineAlignmentCostTable,
};

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchLowerBoundMatrix {
    matrix: Array2<Cost>,
    shift_x: isize,
    shift_y: isize,
    min_distance_between_two_template_switches: usize,
}

type TSLBAlignmentStrategies<AlphabetType> = AlignmentStrategySelection<
    AlphabetType,
    CostOnlyNodeOrdStrategy,
    NoTemplateSwitchMinLengthStrategy,
    NoChainingStrategy,
    MaxTemplateSwitchCountStrategy,
    ForbidSecondaryDeletionStrategy,
    NoShortcutStrategy,
>;

impl TemplateSwitchLowerBoundMatrix {
    pub fn new<AlphabetType: Alphabet>(config: &TemplateSwitchConfig<AlphabetType>) -> Self {
        let lower_bound_config = generate_template_switch_lower_bound_config(config);
        assert!(
            lower_bound_config
                .secondary_edit_costs
                .min_gap_extend_cost()
                > Cost::ZERO,
            "Secondary gap extend costs must be greater than zero for all alphabet characters."
        );

        let mut open_lower_bounds = HashSet::new();
        open_lower_bounds.insert((0isize, 0isize));
        let mut closed_lower_bounds = HashMap::new();
        let mut previous_closed_lower_bounds = HashMap::new();
        let mut genome_length = 1_000;

        'outer: loop {
            debug!("Using genome length {genome_length}");
            assert!(genome_length < usize::try_from(isize::MAX).unwrap() / 2);
            let genome = VectorGenome::<AlphabetType>::from_iter(
                iter::repeat(AlphabetType::iter().next().unwrap()).take(genome_length),
            );
            let mut a_star = AStar::new(Context::<_, TSLBAlignmentStrategies<AlphabetType>>::new(
                genome.as_genome_subsequence(),
                genome.as_genome_subsequence(),
                lower_bound_config.clone(),
                1,
            ));
            let root_xy = genome_length / 2;
            a_star.initialise_with(|context| Node::new_root_at(root_xy, root_xy, context));
            previous_closed_lower_bounds.extend(closed_lower_bounds.drain());

            let root_xy_isize = isize::try_from(root_xy).unwrap();
            let string_length_isize = isize::try_from(genome_length).unwrap();

            'inner: while let Some(coordinates) = open_lower_bounds.iter().next().copied() {
                let (x, y) = open_lower_bounds.take(&coordinates).unwrap();
                if closed_lower_bounds.contains_key(&(x, y)) {
                    continue 'inner;
                }

                trace!("Searching for target ({x}, {y})");

                let has_target = match a_star.search_until(|context, node| {
                    match *node.identifier() {
                        Identifier::Primary {
                            reference_index,
                            query_index,
                            ..
                        } => {
                            reference_index == 0
                                || query_index == 0
                                || reference_index == genome_length - 1
                                || query_index == genome_length - 1
                        }
                        Identifier::PrimaryReentry {
                            reference_index,
                            query_index,
                            ..
                        } => {
                            let reentry_x =
                                isize::try_from(reference_index).unwrap() - root_xy_isize;
                            let reentry_y = isize::try_from(query_index).unwrap() - root_xy_isize;

                            trace!("Found target {} at cost {}", node.identifier(), node.cost());

                            let previous =
                                closed_lower_bounds.insert((reentry_x, reentry_y), node.cost());

                            // If previous exists, then it was inserted when searching with a shorter string.
                            debug_assert!(previous.is_none());
                            let previous_previous = previous_closed_lower_bounds
                                .get(&(reentry_x, reentry_y))
                                .copied();
                            debug_assert!(
                                previous_previous.is_none()
                                    || previous_previous == Some(node.cost())
                            );

                            enqueue_neighbours(
                                reentry_x,
                                reentry_y,
                                &mut closed_lower_bounds,
                                &mut open_lower_bounds,
                            );

                            /*trace!(
                                "Found primary reentry node {} at cost {}",
                                node.identifier(),
                                node.cost()
                            );*/
                            reentry_x == x && reentry_y == y
                        }
                        Identifier::TemplateSwitchEntrance {
                            entrance_reference_index,
                            entrance_query_index,
                            template_switch_first_offset,
                            ..
                        } => {
                            entrance_reference_index == 0
                                || entrance_query_index == 0
                                || entrance_reference_index == genome_length - 1
                                || entrance_query_index == genome_length - 1
                                || template_switch_first_offset + root_xy_isize == 0
                                || template_switch_first_offset + root_xy_isize
                                    == string_length_isize - 1
                        }
                        Identifier::Secondary {
                            primary_index,
                            secondary_index,
                            ..
                        } => {
                            primary_index == 0
                                || primary_index == genome_length - 1
                                || secondary_index == 0
                                || secondary_index == genome_length - 1
                        }
                        Identifier::TemplateSwitchExit { .. } => {
                            node.generate_primary_reentry_successor(context).is_none()
                        }
                    }
                }) {
                    AStarResult::FoundTarget { identifier, cost } => {
                        trace!("Search termianted with target {identifier} at cost {cost}");

                        if let Identifier::PrimaryReentry { .. } = identifier {
                            debug_assert_eq!(closed_lower_bounds.get(&(x, y)), Some(cost).as_ref());
                            true
                        } else {
                            trace!("{:?}", {
                                let mut alignment = Vec::new();

                                // Backtrack.
                                for alignment_type in a_star.backtrack().into_iter().map(
                                    <Context<
                                        SliceSubGenome<_>,
                                        TSLBAlignmentStrategies<AlphabetType>,
                                    > as AlignmentContext>::AlignmentType::from,
                                ) {
                                    if let Some((count, previous_alignment_type)) =
                                        alignment.last_mut()
                                    {
                                        if alignment_type.is_repeated(previous_alignment_type) {
                                            *count += 1;
                                        } else {
                                            alignment.push((1, alignment_type));
                                        }
                                    } else {
                                        alignment.push((1, alignment_type));
                                    }
                                }

                                alignment.reverse();
                                alignment
                            });
                            open_lower_bounds.insert((x, y));
                            genome_length *= 2;
                            continue 'outer;
                        }
                    }
                    AStarResult::NoTarget => {
                        trace!("Search terminated without target");
                        let previous = closed_lower_bounds.insert((x, y), Cost::MAX);
                        debug_assert!(previous.is_none());
                        false
                    }
                };

                if has_target {
                    trace!("Inserting neighbours after search");
                    enqueue_neighbours(x, y, &mut closed_lower_bounds, &mut open_lower_bounds);
                }
            }

            break;
        }

        let min_distance_between_two_template_switches =
            usize::try_from(config.left_flank_length + config.right_flank_length).unwrap();

        let (min_x, max_x, min_y, max_y) = closed_lower_bounds.iter().fold(
            (isize::MAX, isize::MIN, isize::MAX, isize::MIN),
            |(min_x, max_x, min_y, max_y), (&(x, y), &cost)| {
                if cost != Cost::MAX {
                    (min_x.min(x), max_x.max(x), min_y.min(y), max_y.max(y))
                } else {
                    (min_x, max_x, min_y, max_y)
                }
            },
        );

        let shift_x = -min_x;
        let shift_y = -min_y;
        let len_x = usize::try_from(max_x - min_x + 1).unwrap();
        let len_y = usize::try_from(max_y - min_y + 1).unwrap();

        let mut matrix = Array2::from_elem([len_x, len_y], Cost::MAX);
        for ((x, y), cost) in closed_lower_bounds {
            if cost != Cost::MAX {
                let x = usize::try_from(x + shift_x).unwrap();
                let y = usize::try_from(y + shift_y).unwrap();
                matrix[(x, y)] = cost;
            }
        }

        Self {
            matrix,
            shift_x,
            shift_y,
            min_distance_between_two_template_switches,
        }
    }

    pub fn min_distance_between_two_template_switches(&self) -> usize {
        self.min_distance_between_two_template_switches
    }
}

fn generate_template_switch_lower_bound_config<AlphabetType: Alphabet>(
    config: &TemplateSwitchConfig<AlphabetType>,
) -> TemplateSwitchConfig<AlphabetType> {
    TemplateSwitchConfig {
        left_flank_length: 0,
        right_flank_length: 0,
        min_length: config.min_length,

        base_cost: config.base_cost,

        primary_edit_costs: GapAffineAlignmentCostTable::new_max(),
        secondary_edit_costs: config.secondary_edit_costs.clone().into_lower_bound(),
        left_flank_edit_costs: GapAffineAlignmentCostTable::new_max(),
        right_flank_edit_costs: GapAffineAlignmentCostTable::new_max(),

        // The offset only affects which part of the secondary string is being compared against, but otherwise does not change anything.
        // Hence we can ignore it for the lower bound, and simply choose its minimum.
        offset_costs: vec![
            (isize::MIN, Cost::MAX),
            (0, config.offset_costs.min(..).unwrap()),
            (1, Cost::MAX),
        ]
        .try_into()
        .unwrap(),
        length_costs: config.length_costs.clone(),
        length_difference_costs: config.length_difference_costs.clone(),
    }
}

fn enqueue_neighbours(
    x: isize,
    y: isize,
    closed_lower_bounds: &mut HashMap<(isize, isize), Cost>,
    open_lower_bounds: &mut HashSet<(isize, isize)>,
) {
    for (x, y) in [
        (x + 1, y + 1),
        (x + 1, y),
        (x + 1, y - 1),
        (x, y + 1),
        (x, y - 1),
        (x - 1, y + 1),
        (x - 1, y),
        (x - 1, y - 1),
    ] {
        if !closed_lower_bounds.contains_key(&(x, y)) {
            open_lower_bounds.insert((x, y));
        }
    }
}

impl Display for TemplateSwitchLowerBoundMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "TemplateSwitchLowerBoundMatrix:")?;
        writeln!(
            f,
            "min distance: {}",
            self.min_distance_between_two_template_switches
        )?;
        writeln!(
            f,
            "x range: [{}, {}]; y range: [{}, {}]",
            -self.shift_x,
            isize::try_from(self.matrix.dim().0).unwrap() - self.shift_x,
            -self.shift_y,
            isize::try_from(self.matrix.dim().1).unwrap() - self.shift_y
        )?;

        let mut buffer = String::new();
        let column_widths: Vec<_> = (0..self.matrix.dim().0)
            .map(|x| {
                (0..self.matrix.dim().1)
                    .map(|y| {
                        if self.matrix[(x, y)] == Cost::MAX {
                            1
                        } else {
                            buffer.clear();
                            write!(buffer, "{}", self.matrix[(x, y)]).unwrap();
                            buffer.len()
                        }
                    })
                    .max()
                    .unwrap()
            })
            .collect();

        for y in 0..self.matrix.dim().1 {
            for (x, column_width) in column_widths.iter().enumerate().take(self.matrix.dim().0) {
                buffer.clear();
                if self.matrix[(x, y)] == Cost::MAX {
                    write!(buffer, "∞").unwrap();
                } else {
                    write!(buffer, "{}", self.matrix[(x, y)]).unwrap();
                }
                for _ in 0..=column_width - buffer.chars().count() {
                    write!(f, " ")?;
                }
                write!(f, "{buffer}")?;
            }

            writeln!(f)?;
        }

        Ok(())
    }
}
