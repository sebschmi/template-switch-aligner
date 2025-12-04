use generic_a_star::cost::AStarCost;
use ndarray::{Array, Array2};

use crate::{
    alignment::ts_kind::TsKind, anchors::Anchors, chaining_lower_bounds::ChainingLowerBounds,
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
    ) -> Self {
        let k = usize::try_from(chaining_lower_bounds.max_match_run() + 1).unwrap();

        let mut primary = Array2::from_elem(
            (anchors.primary.len(), anchors.primary.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(to_anchor, k) {
                    primary[[from_index, to_index]] =
                        chaining_lower_bounds.primary_lower_bound(gap1, gap2);
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
            }
        }

        let mut jump_12_to_11 = Array::from_elem(
            (anchors.primary.len(), anchors.secondary_11.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_11.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS11, k) {
                    jump_12_to_11[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_12_to_12 = Array::from_elem(
            (anchors.primary.len(), anchors.secondary_12.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_12.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS12, k) {
                    jump_12_to_12[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_12_to_21 = Array::from_elem(
            (anchors.primary.len(), anchors.secondary_21.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_21.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS21, k) {
                    jump_12_to_21[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_12_to_22 = Array::from_elem(
            (anchors.primary.len(), anchors.secondary_22.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.primary.iter().enumerate() {
            for (to_index, to_anchor) in anchors.secondary_22.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS22, k) {
                    jump_12_to_22[[from_index, to_index]] =
                        chaining_lower_bounds.jump_12_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_11 = Array::from_elem(
            (anchors.secondary_11.len(), anchors.primary.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_11.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS11, k) {
                    jump_34_from_11[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_12 = Array::from_elem(
            (anchors.secondary_12.len(), anchors.primary.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_12.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS12, k) {
                    jump_34_from_12[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_21 = Array::from_elem(
            (anchors.secondary_21.len(), anchors.primary.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_21.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS21, k) {
                    jump_34_from_21[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
        }

        let mut jump_34_from_22 = Array::from_elem(
            (anchors.secondary_22.len(), anchors.primary.len()),
            Cost::max_value(),
        );
        for (from_index, from_anchor) in anchors.secondary_22.iter().enumerate() {
            for (to_index, to_anchor) in anchors.primary.iter().enumerate() {
                if let Some(gap) = from_anchor.chaining_jump_gap(to_anchor, TsKind::TS22, k) {
                    jump_34_from_22[[from_index, to_index]] =
                        chaining_lower_bounds.jump_34_lower_bound(gap);
                }
            }
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
}
