use generic_a_star::cost::AStarCost;

use crate::{
    alignment::{
        coordinates::AlignmentCoordinates,
        ts_kind::{TsAncestor, TsDescendant, TsKind},
    },
    anchors::Anchors,
    chaining_cost_function::cost_array::ChainingCostArray,
    chaining_lower_bounds::ChainingLowerBounds,
};

mod cost_array;

pub struct ChainingCostFunction<Cost> {
    primary: ChainingCostArray<Cost>,
    secondary_11: ChainingCostArray<Cost>,
    secondary_12: ChainingCostArray<Cost>,
    secondary_21: ChainingCostArray<Cost>,
    secondary_22: ChainingCostArray<Cost>,
    jump_12_to_11: ChainingCostArray<Cost>,
    jump_12_to_12: ChainingCostArray<Cost>,
    jump_12_to_21: ChainingCostArray<Cost>,
    jump_12_to_22: ChainingCostArray<Cost>,
    jump_34_from_11: ChainingCostArray<Cost>,
    jump_34_from_12: ChainingCostArray<Cost>,
    jump_34_from_21: ChainingCostArray<Cost>,
    jump_34_from_22: ChainingCostArray<Cost>,
}

impl<Cost: AStarCost> ChainingCostFunction<Cost> {
    pub fn new_from_lower_bounds(
        chaining_lower_bounds: &ChainingLowerBounds<Cost>,
        anchors: &Anchors,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> Self {
        let k = usize::try_from(chaining_lower_bounds.max_match_run() + 1).unwrap();

        let mut primary = ChainingCostArray::new_from_cost(
            [anchors.primary.len() + 2, anchors.primary.len() + 2],
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

        let mut secondary_11 = ChainingCostArray::new_from_cost(
            [anchors.secondary_11.len(), anchors.secondary_11.len()],
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

        let mut secondary_12 = ChainingCostArray::new_from_cost(
            [anchors.secondary_12.len(), anchors.secondary_12.len()],
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

        let mut secondary_21 = ChainingCostArray::new_from_cost(
            [anchors.secondary_21.len(), anchors.secondary_21.len()],
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

        let mut secondary_22 = ChainingCostArray::new_from_cost(
            [anchors.secondary_22.len(), anchors.secondary_22.len()],
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

        let mut jump_12_to_11 = ChainingCostArray::new_from_cost(
            [anchors.primary.len() + 2, anchors.secondary_11.len()],
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

        let mut jump_12_to_12 = ChainingCostArray::new_from_cost(
            [anchors.primary.len() + 2, anchors.secondary_12.len()],
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

        let mut jump_12_to_21 = ChainingCostArray::new_from_cost(
            [anchors.primary.len() + 2, anchors.secondary_21.len()],
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

        let mut jump_12_to_22 = ChainingCostArray::new_from_cost(
            [anchors.primary.len() + 2, anchors.secondary_22.len()],
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

        let mut jump_34_from_11 = ChainingCostArray::new_from_cost(
            [anchors.secondary_11.len(), anchors.primary.len() + 2],
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

        let mut jump_34_from_12 = ChainingCostArray::new_from_cost(
            [anchors.secondary_12.len(), anchors.primary.len() + 2],
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

        let mut jump_34_from_21 = ChainingCostArray::new_from_cost(
            [anchors.secondary_21.len(), anchors.primary.len() + 2],
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

        let mut jump_34_from_22 = ChainingCostArray::new_from_cost(
            [anchors.secondary_22.len(), anchors.primary.len() + 2],
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

        for (index, anchor) in anchors.secondary_11.iter().enumerate() {
            let gap = anchor.chaining_jump_gap_from_start(start, TsKind::TS11);
            jump_12_to_11[[0, index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = anchor.chaining_jump_gap_to_end(end, TsKind::TS11, k);
            jump_34_from_11[[index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.jump_34_lower_bound(gap);
        }

        for (index, anchor) in anchors.secondary_12.iter().enumerate() {
            let gap = anchor.chaining_jump_gap_from_start(start, TsKind::TS12);
            jump_12_to_12[[0, index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = anchor.chaining_jump_gap_to_end(end, TsKind::TS12, k);
            jump_34_from_12[[index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.jump_34_lower_bound(gap);
        }

        for (index, anchor) in anchors.secondary_21.iter().enumerate() {
            let gap = anchor.chaining_jump_gap_from_start(start, TsKind::TS21);
            jump_12_to_21[[0, index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = anchor.chaining_jump_gap_to_end(end, TsKind::TS21, k);
            jump_34_from_21[[index, anchors.primary.len() + 1]] =
                chaining_lower_bounds.jump_34_lower_bound(gap);
        }

        for (index, anchor) in anchors.secondary_22.iter().enumerate() {
            let gap = anchor.chaining_jump_gap_from_start(start, TsKind::TS22);
            jump_12_to_22[[0, index]] = chaining_lower_bounds.jump_12_lower_bound(gap);
            let gap = anchor.chaining_jump_gap_to_end(end, TsKind::TS22, k);
            jump_34_from_22[[index, anchors.primary.len() + 1]] =
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

    pub fn jump_12_to_11(&self, from_primary_index: usize, to_secondary_11_index: usize) -> Cost {
        self.jump_12_to_11[[from_primary_index + 1, to_secondary_11_index]]
    }

    pub fn jump_12_to_12(&self, from_primary_index: usize, to_secondary_12_index: usize) -> Cost {
        self.jump_12_to_12[[from_primary_index + 1, to_secondary_12_index]]
    }

    pub fn jump_12_to_21(&self, from_primary_index: usize, to_secondary_21_index: usize) -> Cost {
        self.jump_12_to_21[[from_primary_index + 1, to_secondary_21_index]]
    }

    pub fn jump_12_to_22(&self, from_primary_index: usize, to_secondary_22_index: usize) -> Cost {
        self.jump_12_to_22[[from_primary_index + 1, to_secondary_22_index]]
    }

    pub fn jump_12(
        &self,
        from_primary_index: usize,
        to_secondary_index: usize,
        ts_kind: TsKind,
    ) -> Cost {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.jump_12_to_11(from_primary_index, to_secondary_index)
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.jump_12_to_12(from_primary_index, to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.jump_12_to_21(from_primary_index, to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.jump_12_to_22(from_primary_index, to_secondary_index)
            }
        }
    }

    pub fn jump_12_to_11_from_start(&self, to_secondary_11_index: usize) -> Cost {
        self.jump_12_to_11[[0, to_secondary_11_index]]
    }

    pub fn jump_12_to_12_from_start(&self, to_secondary_12_index: usize) -> Cost {
        self.jump_12_to_12[[0, to_secondary_12_index]]
    }

    pub fn jump_12_to_21_from_start(&self, to_secondary_21_index: usize) -> Cost {
        self.jump_12_to_21[[0, to_secondary_21_index]]
    }

    pub fn jump_12_to_22_from_start(&self, to_secondary_22_index: usize) -> Cost {
        self.jump_12_to_22[[0, to_secondary_22_index]]
    }

    pub fn jump_12_from_start(&self, to_secondary_index: usize, ts_kind: TsKind) -> Cost {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.jump_12_to_11_from_start(to_secondary_index)
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.jump_12_to_12_from_start(to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.jump_12_to_21_from_start(to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.jump_12_to_22_from_start(to_secondary_index)
            }
        }
    }

    pub fn secondary(
        &self,
        from_secondary_index: usize,
        to_secondary_index: usize,
        ts_kind: TsKind,
    ) -> Cost {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.secondary_11[[from_secondary_index, to_secondary_index]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.secondary_12[[from_secondary_index, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.secondary_21[[from_secondary_index, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.secondary_22[[from_secondary_index, to_secondary_index]]
            }
        }
    }

    pub fn jump_34(
        &self,
        from_secondary_index: usize,
        to_primary_index: usize,
        ts_kind: TsKind,
    ) -> Cost {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.jump_34_from_11[[from_secondary_index, to_primary_index + 1]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.jump_34_from_12[[from_secondary_index, to_primary_index + 1]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.jump_34_from_21[[from_secondary_index, to_primary_index + 1]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.jump_34_from_22[[from_secondary_index, to_primary_index + 1]]
            }
        }
    }

    pub fn jump_34_to_end(&self, from_secondary_index: usize, ts_kind: TsKind) -> Cost {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.jump_34_from_11[[from_secondary_index, self.jump_34_from_11.dim().1 - 1]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.jump_34_from_12[[from_secondary_index, self.jump_34_from_12.dim().1 - 1]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.jump_34_from_21[[from_secondary_index, self.jump_34_from_21.dim().1 - 1]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.jump_34_from_22[[from_secondary_index, self.jump_34_from_22.dim().1 - 1]]
            }
        }
    }

    pub fn is_primary_exact(&self, from_primary_index: usize, to_primary_index: usize) -> bool {
        self.primary
            .is_exact(from_primary_index + 1, to_primary_index + 1)
    }

    pub fn is_primary_from_start_exact(&self, primary_index: usize) -> bool {
        self.primary.is_exact(0, primary_index + 1)
    }

    pub fn is_primary_to_end_exact(&self, primary_index: usize) -> bool {
        self.primary
            .is_exact(primary_index + 1, self.primary.dim().1 - 1)
    }

    pub fn is_start_to_end_exact(&self) -> bool {
        self.primary.is_exact(0, self.primary.dim().1 - 1)
    }

    pub fn is_jump_12_to_11_exact(
        &self,
        from_primary_index: usize,
        to_secondary_11_index: usize,
    ) -> bool {
        self.jump_12_to_11
            .is_exact(from_primary_index + 1, to_secondary_11_index)
    }

    pub fn is_jump_12_to_12_exact(
        &self,
        from_primary_index: usize,
        to_secondary_12_index: usize,
    ) -> bool {
        self.jump_12_to_12
            .is_exact(from_primary_index + 1, to_secondary_12_index)
    }

    pub fn is_jump_12_to_21_exact(
        &self,
        from_primary_index: usize,
        to_secondary_21_index: usize,
    ) -> bool {
        self.jump_12_to_21
            .is_exact(from_primary_index + 1, to_secondary_21_index)
    }

    pub fn is_jump_12_to_22_exact(
        &self,
        from_primary_index: usize,
        to_secondary_22_index: usize,
    ) -> bool {
        self.jump_12_to_22
            .is_exact(from_primary_index + 1, to_secondary_22_index)
    }

    pub fn is_jump_12_exact(
        &self,
        from_primary_index: usize,
        to_secondary_index: usize,
        ts_kind: TsKind,
    ) -> bool {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.is_jump_12_to_11_exact(from_primary_index, to_secondary_index)
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.is_jump_12_to_12_exact(from_primary_index, to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.is_jump_12_to_21_exact(from_primary_index, to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.is_jump_12_to_22_exact(from_primary_index, to_secondary_index)
            }
        }
    }

    pub fn is_jump_12_to_11_from_start_exact(&self, to_secondary_11_index: usize) -> bool {
        self.jump_12_to_11.is_exact(0, to_secondary_11_index)
    }

    pub fn is_jump_12_to_12_from_start_exact(&self, to_secondary_12_index: usize) -> bool {
        self.jump_12_to_12.is_exact(0, to_secondary_12_index)
    }

    pub fn is_jump_12_to_21_from_start_exact(&self, to_secondary_21_index: usize) -> bool {
        self.jump_12_to_21.is_exact(0, to_secondary_21_index)
    }

    pub fn is_jump_12_to_22_from_start_exact(&self, to_secondary_22_index: usize) -> bool {
        self.jump_12_to_22.is_exact(0, to_secondary_22_index)
    }

    pub fn is_jump_12_from_start_exact(&self, to_secondary_index: usize, ts_kind: TsKind) -> bool {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.is_jump_12_to_11_from_start_exact(to_secondary_index)
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.is_jump_12_to_12_from_start_exact(to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.is_jump_12_to_21_from_start_exact(to_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.is_jump_12_to_22_from_start_exact(to_secondary_index)
            }
        }
    }

    pub fn is_secondary_exact(
        &self,
        from_secondary_index: usize,
        to_secondary_index: usize,
        ts_kind: TsKind,
    ) -> bool {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => self
                .secondary_11
                .is_exact(from_secondary_index, to_secondary_index),
            (TsAncestor::Seq1, TsDescendant::Seq2) => self
                .secondary_12
                .is_exact(from_secondary_index, to_secondary_index),
            (TsAncestor::Seq2, TsDescendant::Seq1) => self
                .secondary_21
                .is_exact(from_secondary_index, to_secondary_index),
            (TsAncestor::Seq2, TsDescendant::Seq2) => self
                .secondary_22
                .is_exact(from_secondary_index, to_secondary_index),
        }
    }

    pub fn is_jump_34_exact(
        &self,
        from_secondary_index: usize,
        to_primary_index: usize,
        ts_kind: TsKind,
    ) -> bool {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => self
                .jump_34_from_11
                .is_exact(from_secondary_index, to_primary_index + 1),
            (TsAncestor::Seq1, TsDescendant::Seq2) => self
                .jump_34_from_12
                .is_exact(from_secondary_index, to_primary_index + 1),
            (TsAncestor::Seq2, TsDescendant::Seq1) => self
                .jump_34_from_21
                .is_exact(from_secondary_index, to_primary_index + 1),
            (TsAncestor::Seq2, TsDescendant::Seq2) => self
                .jump_34_from_22
                .is_exact(from_secondary_index, to_primary_index + 1),
        }
    }

    pub fn is_jump_34_to_end_exact(&self, from_secondary_index: usize, ts_kind: TsKind) -> bool {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => self
                .jump_34_from_11
                .is_exact(from_secondary_index, self.jump_34_from_11.dim().1 - 1),
            (TsAncestor::Seq1, TsDescendant::Seq2) => self
                .jump_34_from_12
                .is_exact(from_secondary_index, self.jump_34_from_12.dim().1 - 1),
            (TsAncestor::Seq2, TsDescendant::Seq1) => self
                .jump_34_from_21
                .is_exact(from_secondary_index, self.jump_34_from_21.dim().1 - 1),
            (TsAncestor::Seq2, TsDescendant::Seq2) => self
                .jump_34_from_22
                .is_exact(from_secondary_index, self.jump_34_from_22.dim().1 - 1),
        }
    }

    pub fn update_primary(
        &mut self,
        from_primary_index: usize,
        to_primary_index: usize,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.primary
                .set_exact(from_primary_index + 1, to_primary_index + 1);
        }
        let target = &mut self.primary[[from_primary_index + 1, to_primary_index + 1]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_primary_from_start(
        &mut self,
        primary_index: usize,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.primary.set_exact(0, primary_index + 1);
        }
        let target = &mut self.primary[[0, primary_index + 1]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_primary_to_end(
        &mut self,
        primary_index: usize,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let end_index = self.primary.dim().1 - 1;
        if is_exact {
            self.primary.set_exact(primary_index + 1, end_index);
        }
        let target = &mut self.primary[[primary_index + 1, end_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_start_to_end(&mut self, cost: Cost, is_exact: bool) -> bool {
        let end_index = self.primary.dim().1 - 1;
        if is_exact {
            self.primary.set_exact(0, end_index);
        }
        let target = &mut self.primary[[0, end_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_12(
        &mut self,
        from_primary_index: usize,
        to_secondary_index: usize,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let target = match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                if is_exact {
                    self.jump_12_to_11
                        .set_exact(from_primary_index + 1, to_secondary_index);
                }
                &mut self.jump_12_to_11[[from_primary_index + 1, to_secondary_index]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                if is_exact {
                    self.jump_12_to_12
                        .set_exact(from_primary_index + 1, to_secondary_index);
                }
                &mut self.jump_12_to_12[[from_primary_index + 1, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                if is_exact {
                    self.jump_12_to_21
                        .set_exact(from_primary_index + 1, to_secondary_index);
                }
                &mut self.jump_12_to_21[[from_primary_index + 1, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                if is_exact {
                    self.jump_12_to_22
                        .set_exact(from_primary_index + 1, to_secondary_index);
                }
                &mut self.jump_12_to_22[[from_primary_index + 1, to_secondary_index]]
            }
        };
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_12_from_start(
        &mut self,
        to_secondary_index: usize,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let target = match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                if is_exact {
                    self.jump_12_to_11.set_exact(0, to_secondary_index);
                }
                &mut self.jump_12_to_11[[0, to_secondary_index]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                if is_exact {
                    self.jump_12_to_12.set_exact(0, to_secondary_index);
                }
                &mut self.jump_12_to_12[[0, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                if is_exact {
                    self.jump_12_to_21.set_exact(0, to_secondary_index);
                }
                &mut self.jump_12_to_21[[0, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                if is_exact {
                    self.jump_12_to_22.set_exact(0, to_secondary_index);
                }
                &mut self.jump_12_to_22[[0, to_secondary_index]]
            }
        };
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_secondary(
        &mut self,
        from_secondary_index: usize,
        to_secondary_index: usize,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let target = match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                if is_exact {
                    self.secondary_11
                        .set_exact(from_secondary_index, to_secondary_index);
                }
                &mut self.secondary_11[[from_secondary_index, to_secondary_index]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                if is_exact {
                    self.secondary_12
                        .set_exact(from_secondary_index, to_secondary_index);
                }
                &mut self.secondary_12[[from_secondary_index, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                if is_exact {
                    self.secondary_21
                        .set_exact(from_secondary_index, to_secondary_index);
                }
                &mut self.secondary_21[[from_secondary_index, to_secondary_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                if is_exact {
                    self.secondary_22
                        .set_exact(from_secondary_index, to_secondary_index);
                }
                &mut self.secondary_22[[from_secondary_index, to_secondary_index]]
            }
        };
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_34(
        &mut self,
        from_secondary_index: usize,
        to_primary_index: usize,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let target = match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                if is_exact {
                    self.jump_34_from_11
                        .set_exact(from_secondary_index, to_primary_index + 1);
                }
                &mut self.jump_34_from_11[[from_secondary_index, to_primary_index + 1]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                if is_exact {
                    self.jump_34_from_12
                        .set_exact(from_secondary_index, to_primary_index + 1);
                }
                &mut self.jump_34_from_12[[from_secondary_index, to_primary_index + 1]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                if is_exact {
                    self.jump_34_from_21
                        .set_exact(from_secondary_index, to_primary_index + 1);
                }
                &mut self.jump_34_from_21[[from_secondary_index, to_primary_index + 1]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                if is_exact {
                    self.jump_34_from_22
                        .set_exact(from_secondary_index, to_primary_index + 1);
                }
                &mut self.jump_34_from_22[[from_secondary_index, to_primary_index + 1]]
            }
        };
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_34_to_end(
        &mut self,
        from_secondary_index: usize,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let target = match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                let end_index = self.jump_34_from_11.dim().1 - 1;
                if is_exact {
                    self.jump_34_from_11
                        .set_exact(from_secondary_index, end_index);
                }
                &mut self.jump_34_from_11[[from_secondary_index, end_index]]
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                let end_index = self.jump_34_from_12.dim().1 - 1;
                if is_exact {
                    self.jump_34_from_12
                        .set_exact(from_secondary_index, end_index);
                }
                &mut self.jump_34_from_12[[from_secondary_index, end_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                let end_index = self.jump_34_from_21.dim().1 - 1;
                if is_exact {
                    self.jump_34_from_21
                        .set_exact(from_secondary_index, end_index);
                }
                &mut self.jump_34_from_21[[from_secondary_index, end_index]]
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                let end_index = self.jump_34_from_22.dim().1 - 1;
                if is_exact {
                    self.jump_34_from_22
                        .set_exact(from_secondary_index, end_index);
                }
                &mut self.jump_34_from_22[[from_secondary_index, end_index]]
            }
        };
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn iter_primary_in_cost_order(
        &mut self,
        from_primary_index: usize,
    ) -> impl Iterator<Item = (usize, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.primary.iter_in_cost_order(from_primary_index + 1)
    }

    pub fn iter_primary_from_start_in_cost_order(&mut self) -> impl Iterator<Item = (usize, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.primary.iter_in_cost_order(0)
    }

    pub fn iter_jump_12_in_cost_order(
        &mut self,
        from_primary_index: usize,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (usize, Cost)>
    where
        Cost: Copy + Ord,
    {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => self
                .jump_12_to_11
                .iter_in_cost_order(from_primary_index + 1),
            (TsAncestor::Seq1, TsDescendant::Seq2) => self
                .jump_12_to_12
                .iter_in_cost_order(from_primary_index + 1),
            (TsAncestor::Seq2, TsDescendant::Seq1) => self
                .jump_12_to_21
                .iter_in_cost_order(from_primary_index + 1),
            (TsAncestor::Seq2, TsDescendant::Seq2) => self
                .jump_12_to_22
                .iter_in_cost_order(from_primary_index + 1),
        }
    }

    pub fn iter_jump_12_from_start_in_cost_order(
        &mut self,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (usize, Cost)>
    where
        Cost: Copy + Ord,
    {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => self.jump_12_to_11.iter_in_cost_order(0),
            (TsAncestor::Seq1, TsDescendant::Seq2) => self.jump_12_to_12.iter_in_cost_order(0),
            (TsAncestor::Seq2, TsDescendant::Seq1) => self.jump_12_to_21.iter_in_cost_order(0),
            (TsAncestor::Seq2, TsDescendant::Seq2) => self.jump_12_to_22.iter_in_cost_order(0),
        }
    }

    pub fn iter_secondary_in_cost_order(
        &mut self,
        from_secondary_index: usize,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (usize, Cost)>
    where
        Cost: Copy + Ord,
    {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => {
                self.secondary_11.iter_in_cost_order(from_secondary_index)
            }
            (TsAncestor::Seq1, TsDescendant::Seq2) => {
                self.secondary_12.iter_in_cost_order(from_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq1) => {
                self.secondary_21.iter_in_cost_order(from_secondary_index)
            }
            (TsAncestor::Seq2, TsDescendant::Seq2) => {
                self.secondary_22.iter_in_cost_order(from_secondary_index)
            }
        }
    }

    pub fn iter_jump_34_in_cost_order(
        &mut self,
        from_secondary_index: usize,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (usize, Cost)>
    where
        Cost: Copy + Ord,
    {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => self
                .jump_34_from_11
                .iter_in_cost_order(from_secondary_index),
            (TsAncestor::Seq1, TsDescendant::Seq2) => self
                .jump_34_from_12
                .iter_in_cost_order(from_secondary_index),
            (TsAncestor::Seq2, TsDescendant::Seq1) => self
                .jump_34_from_21
                .iter_in_cost_order(from_secondary_index),
            (TsAncestor::Seq2, TsDescendant::Seq2) => self
                .jump_34_from_22
                .iter_in_cost_order(from_secondary_index),
        }
    }
}
