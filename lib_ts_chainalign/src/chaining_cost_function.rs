use generic_a_star::cost::AStarCost;
use ndarray::{Array, Array2};

use crate::{
    alignment::{coordinates::AlignmentCoordinates, ts_kind::TsKind},
    anchors::Anchors,
    chaining_lower_bounds::ChainingLowerBounds,
};

pub struct ChainingCostFunction<Cost> {
    primary: Array2<Cost>,
    secondary_11: Array2<Cost>,
    secondary_12: Array2<Cost>,
    secondary_21: Array2<Cost>,
    secondary_22: Array2<Cost>,
    jump_12_to_11: Array2<Cost>,
    jump_12_to_12: Array2<Cost>,
    jump_12_to_21: Array2<Cost>,
    jump_12_to_22: Array2<Cost>,
    jump_34_from_11: Array2<Cost>,
    jump_34_from_12: Array2<Cost>,
    jump_34_from_21: Array2<Cost>,
    jump_34_from_22: Array2<Cost>,
}

impl<Cost: AStarCost> ChainingCostFunction<Cost> {
    pub fn new_from_lower_bounds(
        chaining_lower_bounds: &ChainingLowerBounds<Cost>,
        anchors: &Anchors,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> Self {
        let k = usize::try_from(chaining_lower_bounds.max_match_run() + 1).unwrap();

        let mut primary = Array2::from_elem(
            (anchors.primary.len() + 2, anchors.primary.len() + 2),
            Cost::max_value(),
        );
        let gap1 = end.primary_ordinate_a().unwrap() - start.primary_ordinate_a().unwrap();
        let gap2 = end.primary_ordinate_b().unwrap() - start.primary_ordinate_b().unwrap();
        primary[[0, anchors.primary.len() + 1]] =
            chaining_lower_bounds.primary_lower_bound(gap1, gap2);
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            let from_index = from_index + 1;
            let (gap1, gap2) = from_anchor.chaining_gaps_from_start(start);
            primary[[0, from_index]] = chaining_lower_bounds.primary_lower_bound(gap1, gap2);
            let (gap1, gap2) = from_anchor.chaining_gaps_to_end(end, k);
            primary[[from_index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.primary_lower_bound(gap1, gap2);

            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                let to_index = to_index + 1;
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(to_anchor, k) {
                    primary[[from_index, to_index]] =
                        chaining_lower_bounds.primary_lower_bound(gap1, gap2);
                }
                if from_anchor.is_direct_predecessor_of(to_anchor) {
                    primary[[from_index, to_index]] = Cost::zero();
                }
            }
        }

        let mut secondary_11 = Array2::from_elem(
            (anchors.secondary_11.len(), anchors.secondary_11.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_11.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_11.iter().enumerate() {
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(to_anchor, TsKind::TS11, k) {
                    secondary_11[[from_index, to_index]] =
                        chaining_lower_bounds.secondary_lower_bound(gap1, gap2);
                }
                if from_anchor.is_direct_predecessor_of(to_anchor) {
                    secondary_11[[from_index, to_index]] = Cost::zero();
                }
            }
        }

        let mut secondary_12 = Array2::from_elem(
            (anchors.secondary_12.len(), anchors.secondary_12.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_12.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_12.iter().enumerate() {
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(to_anchor, TsKind::TS12, k) {
                    secondary_12[[from_index, to_index]] =
                        chaining_lower_bounds.secondary_lower_bound(gap1, gap2);
                }
                if from_anchor.is_direct_predecessor_of(to_anchor) {
                    secondary_12[[from_index, to_index]] = Cost::zero();
                }
            }
        }

        let mut secondary_21 = Array2::from_elem(
            (anchors.secondary_21.len(), anchors.secondary_21.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_21.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_21.iter().enumerate() {
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(to_anchor, TsKind::TS21, k) {
                    secondary_21[[from_index, to_index]] =
                        chaining_lower_bounds.secondary_lower_bound(gap1, gap2);
                }
                if from_anchor.is_direct_predecessor_of(to_anchor) {
                    secondary_21[[from_index, to_index]] = Cost::zero();
                }
            }
        }

        let mut secondary_22 = Array2::from_elem(
            (anchors.secondary_22.len(), anchors.secondary_22.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_22.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_22.iter().enumerate() {
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(to_anchor, TsKind::TS22, k) {
                    secondary_22[[from_index, to_index]] =
                        chaining_lower_bounds.secondary_lower_bound(gap1, gap2);
                }
                if from_anchor.is_direct_predecessor_of(to_anchor) {
                    secondary_22[[from_index, to_index]] = Cost::zero();
                }
            }
        }

        let mut jump_12_to_11 = Array::from_elem(
            (anchors.primary.len() + 2, anchors.secondary_11.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            let from_index = from_index + 1;
            for (to_index, to_anchor) in anchors.secondary_11.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS11, k) {
                    jump_12_to_11[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_12_to_12 = Array::from_elem(
            (anchors.primary.len() + 2, anchors.secondary_12.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            let from_index = from_index + 1;
            for (to_index, to_anchor) in anchors.secondary_12.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS12, k) {
                    jump_12_to_12[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_12_to_21 = Array::from_elem(
            (anchors.primary.len() + 2, anchors.secondary_21.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            let from_index = from_index + 1;
            for (to_index, to_anchor) in anchors.secondary_21.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS21, k) {
                    jump_12_to_21[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_12_to_22 = Array::from_elem(
            (anchors.primary.len() + 2, anchors.secondary_22.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            let from_index = from_index + 1;
            for (to_index, to_anchor) in anchors.secondary_22.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS22, k) {
                    jump_12_to_22[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_11 = Array::from_elem(
            (anchors.secondary_11.len(), anchors.primary.len() + 2),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_11.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                let to_index = to_index + 1;
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS11, k) {
                    jump_34_from_11[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_12 = Array::from_elem(
            (anchors.secondary_12.len(), anchors.primary.len() + 2),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_12.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                let to_index = to_index + 1;
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS12, k) {
                    jump_34_from_12[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_21 = Array::from_elem(
            (anchors.secondary_21.len(), anchors.primary.len() + 2),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_21.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                let to_index = to_index + 1;
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS21, k) {
                    jump_34_from_21[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_22 = Array::from_elem(
            (anchors.secondary_22.len(), anchors.primary.len() + 2),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_22.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                let to_index = to_index + 1;
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS22, k) {
                    jump_34_from_22[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
        }

        for (to_index, to_anchor) in anchors.secondary_11.iter().enumerate() {
            let gap = to_anchor.chaining_jump_gap_from_start(start, TsKind::TS11);
            jump_12_to_11[[0, to_index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = to_anchor.chaining_jump_gap_to_end(end, TsKind::TS11, k);
            jump_34_from_11[[to_index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.jump_34_lower_bound(gap);
        }

        for (to_index, to_anchor) in anchors.secondary_12.iter().enumerate() {
            let gap = to_anchor.chaining_jump_gap_from_start(start, TsKind::TS12);
            jump_12_to_12[[0, to_index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = to_anchor.chaining_jump_gap_to_end(end, TsKind::TS12, k);
            jump_34_from_12[[to_index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.jump_34_lower_bound(gap);
        }

        for (to_index, to_anchor) in anchors.secondary_21.iter().enumerate() {
            let gap = to_anchor.chaining_jump_gap_from_start(start, TsKind::TS21);
            jump_12_to_21[[0, to_index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = to_anchor.chaining_jump_gap_to_end(end, TsKind::TS21, k);
            jump_34_from_21[[to_index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.jump_34_lower_bound(gap);
        }

        for (to_index, to_anchor) in anchors.secondary_22.iter().enumerate() {
            let gap = to_anchor.chaining_jump_gap_from_start(start, TsKind::TS22);
            jump_12_to_22[[0, to_index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = to_anchor.chaining_jump_gap_to_end(end, TsKind::TS22, k);
            jump_34_from_22[[to_index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.jump_34_lower_bound(gap);
        }

        Self {
            primary,
            secondary_11,
            secondary_12,
            secondary_21,
            secondary_22,
            jump_12_to_11,
            jump_12_to_12,
            jump_12_to_21,
            jump_12_to_22,
            jump_34_from_11,
            jump_34_from_12,
            jump_34_from_21,
            jump_34_from_22,
        }
    }

    pub fn primary(&self, from_primary_index: usize, to_primary_index: usize) -> Cost {
        self.primary[[from_primary_index + 1, to_primary_index + 1]]
    }

    pub fn primary_from_start(&self, primary_index: usize) -> Cost {
        self.primary[[0, primary_index + 1]]
    }

    pub fn primary_to_end(&self, primary_index: usize) -> Cost {
        self.primary[[primary_index + 1, self.primary.dim().1 - 1]]
    }

    pub fn start_to_end(&self) -> Cost {
        self.primary[[0, self.primary.dim().1 - 1]]
    }
}
