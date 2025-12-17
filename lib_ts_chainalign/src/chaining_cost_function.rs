use std::time::Instant;

use generic_a_star::cost::AStarCost;
use itertools::Itertools;
use log::debug;
use num_traits::Zero;

use crate::{
    alignment::{coordinates::AlignmentCoordinates, ts_kind::TsKind},
    anchors::{Anchors, index::AnchorIndex, primary::PrimaryAnchor, secondary::SecondaryAnchor},
    chaining_cost_function::cost_array::ChainingCostArray,
    chaining_lower_bounds::ChainingLowerBounds,
};

mod cost_array;

pub struct ChainingCostFunction<Cost> {
    primary: ChainingCostArray<Cost>,
    secondaries: [ChainingCostArray<Cost>; 4],
    jump_12s: [ChainingCostArray<Cost>; 4],
    jump_34s: [ChainingCostArray<Cost>; 4],
}

impl<Cost: AStarCost> ChainingCostFunction<Cost> {
    pub fn new_from_lower_bounds(
        chaining_lower_bounds: &ChainingLowerBounds<Cost>,
        anchors: &Anchors,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> Self {
        let start_time = Instant::now();

        let k = usize::try_from(chaining_lower_bounds.max_match_run() + 1).unwrap();
        let primary_anchor_amount = anchors.primary_len() + 2;
        let primary_start_anchor_index = AnchorIndex::zero();
        let primary_end_anchor_index = primary_anchor_amount - 1;

        let mut primary = ChainingCostArray::new_from_cost(
            [primary_anchor_amount, primary_anchor_amount],
            Cost::max_value(),
        );
        let gap1 = end.primary_ordinate_a().unwrap() - start.primary_ordinate_a().unwrap();
        let gap2 = end.primary_ordinate_b().unwrap() - start.primary_ordinate_b().unwrap();
        primary[[primary_start_anchor_index, primary_end_anchor_index]] =
            chaining_lower_bounds.primary_lower_bound(gap1, gap2);
        for (from_index, from_anchor) in anchors.enumerate_primaries() {
            let from_index = from_index + 1;
            let (gap1, gap2) = from_anchor.chaining_gaps_from_start(start);
            primary[[primary_start_anchor_index, from_index]] =
                chaining_lower_bounds.primary_lower_bound(gap1, gap2);
            let (gap1, gap2) = from_anchor.chaining_gaps_to_end(end, k);
            primary[[from_index, primary_end_anchor_index]] =
                chaining_lower_bounds.primary_lower_bound(gap1, gap2);

            for (to_index, to_anchor) in anchors.enumerate_primaries() {
                let to_index = to_index + 1;
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(&to_anchor, k) {
                    primary[[from_index, to_index]] =
                        chaining_lower_bounds.primary_lower_bound(gap1, gap2);
                }
                if from_anchor.is_direct_predecessor_of(&to_anchor) {
                    primary[[from_index, to_index]] = Cost::zero();
                }
            }
        }

        let mut secondaries = TsKind::iter()
            .map(|ts_kind| {
                ChainingCostArray::new_from_cost(
                    [
                        anchors.secondary_len(ts_kind),
                        anchors.secondary_len(ts_kind),
                    ],
                    Cost::max_value(),
                )
            })
            .collect_array()
            .unwrap();
        for (ts_kind, secondary) in TsKind::iter().zip(&mut secondaries) {
            for (from_index, from_anchor) in anchors.enumerate_secondaries(ts_kind) {
                for (to_index, to_anchor) in anchors.enumerate_secondaries(ts_kind) {
                    if let Some((gap1, gap2)) = from_anchor.chaining_gaps(&to_anchor, ts_kind, k) {
                        secondary[[from_index, to_index]] =
                            chaining_lower_bounds.secondary_lower_bound(gap1, gap2);
                    }
                    if from_anchor.is_direct_predecessor_of(&to_anchor) {
                        secondary[[from_index, to_index]] = Cost::zero();
                    }
                }
            }
        }

        let mut jump_12s = TsKind::iter()
            .map(|ts_kind| {
                ChainingCostArray::new_from_cost(
                    [primary_anchor_amount, anchors.secondary_len(ts_kind)],
                    Cost::max_value(),
                )
            })
            .collect_array()
            .unwrap();
        for (ts_kind, jump_12) in TsKind::iter().zip(&mut jump_12s) {
            for (from_index, from_anchor) in anchors.enumerate_primaries() {
                let from_index = from_index + 1;
                for (to_index, to_anchor) in anchors.enumerate_secondaries(ts_kind) {
                    if let Some(gap) = from_anchor.chaining_jump_gap(&to_anchor, ts_kind, k) {
                        jump_12[[from_index, to_index]] =
                            chaining_lower_bounds.jump_12_lower_bound(gap);
                    }
                }
            }
        }

        let mut jump_34s = TsKind::iter()
            .map(|ts_kind| {
                ChainingCostArray::new_from_cost(
                    [anchors.secondary_len(ts_kind), primary_anchor_amount],
                    Cost::max_value(),
                )
            })
            .collect_array()
            .unwrap();
        for (ts_kind, jump_34) in TsKind::iter().zip(&mut jump_34s) {
            for (from_index, from_anchor) in anchors.enumerate_secondaries(ts_kind) {
                for (to_index, to_anchor) in anchors.enumerate_primaries() {
                    let to_index = to_index + 1;
                    if let Some(gap) = from_anchor.chaining_jump_gap(&to_anchor, ts_kind, k) {
                        jump_34[[from_index, to_index]] =
                            chaining_lower_bounds.jump_34_lower_bound(gap);
                    }
                }
            }
        }

        for (ts_kind, (jump_12, jump_34)) in
            TsKind::iter().zip(jump_12s.iter_mut().zip(&mut jump_34s))
        {
            for (index, anchor) in anchors.enumerate_secondaries(ts_kind) {
                let gap = anchor.chaining_jump_gap_from_start(start, ts_kind);
                jump_12[[primary_start_anchor_index, index]] =
                    chaining_lower_bounds.jump_12_lower_bound(gap);
                let gap = anchor.chaining_jump_gap_to_end(end, ts_kind, k);
                jump_34[[index, primary_end_anchor_index]] =
                    chaining_lower_bounds.jump_34_lower_bound(gap);
            }
        }

        let end_time = Instant::now();
        let duration = end_time - start_time;
        debug!(
            "Initialising chaining cost function took {:.0}ms",
            duration.as_secs_f64() * 1000.0
        );

        Self {
            primary,
            secondaries,
            jump_12s,
            jump_34s,
        }
    }

    fn primary_start_anchor_index() -> AnchorIndex {
        AnchorIndex::zero()
    }

    fn primary_end_anchor_index(&self) -> AnchorIndex {
        self.primary.dim().0 - 1
    }

    fn secondary_array(&self, ts_kind: TsKind) -> &ChainingCostArray<Cost> {
        &self.secondaries[ts_kind.index()]
    }

    fn secondary_array_mut(&mut self, ts_kind: TsKind) -> &mut ChainingCostArray<Cost> {
        &mut self.secondaries[ts_kind.index()]
    }

    fn jump_12_array(&self, ts_kind: TsKind) -> &ChainingCostArray<Cost> {
        &self.jump_12s[ts_kind.index()]
    }

    fn jump_12_array_mut(&mut self, ts_kind: TsKind) -> &mut ChainingCostArray<Cost> {
        &mut self.jump_12s[ts_kind.index()]
    }

    fn jump_34_array(&self, ts_kind: TsKind) -> &ChainingCostArray<Cost> {
        &self.jump_34s[ts_kind.index()]
    }

    fn jump_34_array_mut(&mut self, ts_kind: TsKind) -> &mut ChainingCostArray<Cost> {
        &mut self.jump_34s[ts_kind.index()]
    }

    pub fn primary(&self, from_primary_index: AnchorIndex, to_primary_index: AnchorIndex) -> Cost {
        self.primary[[from_primary_index + 1, to_primary_index + 1]]
    }

    pub fn primary_from_start(&self, primary_index: AnchorIndex) -> Cost {
        self.primary[[Self::primary_start_anchor_index(), primary_index + 1]]
    }

    pub fn primary_to_end(&self, primary_index: AnchorIndex) -> Cost {
        self.primary[[primary_index + 1, self.primary_end_anchor_index()]]
    }

    pub fn start_to_end(&self) -> Cost {
        self.primary[[
            Self::primary_start_anchor_index(),
            self.primary_end_anchor_index(),
        ]]
    }

    pub fn jump_12(
        &self,
        from_primary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> Cost {
        self.jump_12_array(ts_kind)[[from_primary_index + 1, to_secondary_index]]
    }

    pub fn jump_12_from_start(&self, to_secondary_index: AnchorIndex, ts_kind: TsKind) -> Cost {
        self.jump_12_array(ts_kind)[[Self::primary_start_anchor_index(), to_secondary_index]]
    }

    pub fn secondary(
        &self,
        from_secondary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> Cost {
        self.secondary_array(ts_kind)[[from_secondary_index, to_secondary_index]]
    }

    pub fn jump_34(
        &self,
        from_secondary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> Cost {
        self.jump_34_array(ts_kind)[[from_secondary_index, to_primary_index + 1]]
    }

    pub fn jump_34_to_end(&self, from_secondary_index: AnchorIndex, ts_kind: TsKind) -> Cost {
        self.jump_34_array(ts_kind)[[from_secondary_index, self.primary_end_anchor_index()]]
    }

    pub fn is_primary_exact(
        &self,
        from_primary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
    ) -> bool {
        self.primary
            .is_exact(from_primary_index + 1, to_primary_index + 1)
    }

    pub fn is_primary_from_start_exact(&self, primary_index: AnchorIndex) -> bool {
        self.primary
            .is_exact(Self::primary_start_anchor_index(), primary_index + 1)
    }

    pub fn is_primary_to_end_exact(&self, primary_index: AnchorIndex) -> bool {
        self.primary
            .is_exact(primary_index + 1, self.primary_end_anchor_index())
    }

    pub fn is_start_to_end_exact(&self) -> bool {
        self.primary.is_exact(
            Self::primary_start_anchor_index(),
            self.primary_end_anchor_index(),
        )
    }

    pub fn is_jump_12_exact(
        &self,
        from_primary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_12_array(ts_kind)
            .is_exact(from_primary_index + 1, to_secondary_index)
    }

    pub fn is_jump_12_from_start_exact(
        &self,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_12_array(ts_kind)
            .is_exact(Self::primary_start_anchor_index(), to_secondary_index)
    }

    pub fn is_secondary_exact(
        &self,
        from_secondary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.secondary_array(ts_kind)
            .is_exact(from_secondary_index, to_secondary_index)
    }

    pub fn is_jump_34_exact(
        &self,
        from_secondary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_34_array(ts_kind)
            .is_exact(from_secondary_index, to_primary_index + 1)
    }

    pub fn is_jump_34_to_end_exact(
        &self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_34_array(ts_kind)
            .is_exact(from_secondary_index, self.primary_end_anchor_index())
    }

    pub fn update_primary(
        &mut self,
        from_primary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.primary
                .set_exact(from_primary_index + 1, to_primary_index + 1);
        }
        let target = &mut self.primary[[from_primary_index + 1, to_primary_index + 1]];
        assert!(
            *target <= cost,
            "Target is larger than cost.\ntarget: {target}; cost: {cost}; from_primary_index: {from_primary_index}; to_primary_index: {to_primary_index}; is_exact: {is_exact}"
        );
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_primary_from_start(
        &mut self,
        primary_index: AnchorIndex,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.primary
                .set_exact(Self::primary_start_anchor_index(), primary_index + 1);
        }
        let target = &mut self.primary[[Self::primary_start_anchor_index(), primary_index + 1]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_primary_to_end(
        &mut self,
        primary_index: AnchorIndex,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let end_index = self.primary_end_anchor_index();
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
            self.primary
                .set_exact(Self::primary_start_anchor_index(), end_index);
        }
        let target = &mut self.primary[[Self::primary_start_anchor_index(), end_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_12(
        &mut self,
        from_primary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.jump_12_array_mut(ts_kind)
                .set_exact(from_primary_index + 1, to_secondary_index);
        }
        let target =
            &mut self.jump_12_array_mut(ts_kind)[[from_primary_index + 1, to_secondary_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_12_from_start(
        &mut self,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.jump_12_array_mut(ts_kind)
                .set_exact(Self::primary_start_anchor_index(), to_secondary_index);
        }
        let target = &mut self.jump_12_array_mut(ts_kind)
            [[Self::primary_start_anchor_index(), to_secondary_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_secondary(
        &mut self,
        from_secondary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.secondary_array_mut(ts_kind)
                .set_exact(from_secondary_index, to_secondary_index);
        }
        let target =
            &mut self.secondary_array_mut(ts_kind)[[from_secondary_index, to_secondary_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_34(
        &mut self,
        from_secondary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.jump_34_array_mut(ts_kind)
                .set_exact(from_secondary_index, to_primary_index + 1);
        }
        let target =
            &mut self.jump_34_array_mut(ts_kind)[[from_secondary_index, to_primary_index + 1]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_34_to_end(
        &mut self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let end_index = self.primary_end_anchor_index();
        if is_exact {
            self.jump_34_array_mut(ts_kind)
                .set_exact(from_secondary_index, end_index);
        }
        let target = &mut self.jump_34_array_mut(ts_kind)[[from_secondary_index, end_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn iter_primary_in_cost_order(
        &mut self,
        from_primary_index: AnchorIndex,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.primary
            .iter_in_cost_order_from(from_primary_index + 1, offset)
            .filter(|(to_primary_index, _)| *to_primary_index != AnchorIndex::zero())
            .map(|(to_primary_index, chaining_cost)| (to_primary_index - 1, chaining_cost))
    }

    pub fn iter_primary_from_start_in_cost_order(
        &mut self,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.primary
            .iter_in_cost_order_from(AnchorIndex::zero(), offset)
            .filter(|(to_primary_index, _)| *to_primary_index != AnchorIndex::zero())
            .map(|(to_primary_index, chaining_cost)| (to_primary_index - 1, chaining_cost))
    }

    pub fn iter_jump_12_in_cost_order(
        &mut self,
        from_primary_index: AnchorIndex,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.jump_12_array_mut(ts_kind)
            .iter_in_cost_order_from(from_primary_index + 1, offset)
    }

    pub fn iter_jump_12_from_start_in_cost_order(
        &mut self,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.jump_12_array_mut(ts_kind)
            .iter_in_cost_order_from(Self::primary_start_anchor_index(), offset)
    }

    pub fn iter_secondary_in_cost_order(
        &mut self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.secondary_array_mut(ts_kind)
            .iter_in_cost_order_from(from_secondary_index, offset)
    }

    pub fn iter_jump_34_in_cost_order(
        &mut self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.jump_34_array_mut(ts_kind)
            .iter_in_cost_order_from(from_secondary_index, offset)
            .filter(|(to_primary_index, _)| *to_primary_index != AnchorIndex::zero())
            .map(|(to_primary_index, chaining_cost)| (to_primary_index - 1, chaining_cost))
    }

    pub fn update_additional_primary_targets(
        &mut self,
        from_primary_index: AnchorIndex,
        additional_targets: &mut [(PrimaryAnchor, Cost)],
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_anchor_index, cost) in anchors
            .primary_anchor_to_index_iter(additional_targets.iter().map(|(anchor, _)| *anchor))
            .zip(additional_targets.iter().map(|(_, cost)| *cost))
        {
            let Some(to_primary_index) = to_anchor_index else {
                continue;
            };

            if self.is_primary_exact(from_primary_index, to_primary_index) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(self.primary(from_primary_index, to_primary_index), cost,);
            } else {
                self.update_primary(from_primary_index, to_primary_index, cost, true);
            }
        }
    }

    pub fn update_additional_primary_targets_from_start(
        &mut self,
        additional_targets: &mut [(PrimaryAnchor, Cost)],
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_anchor_index, cost) in anchors
            .primary_anchor_to_index_iter(additional_targets.iter().map(|(anchor, _)| *anchor))
            .zip(additional_targets.iter().map(|(_, cost)| *cost))
        {
            let Some(to_primary_index) = to_anchor_index else {
                continue;
            };

            if self.is_primary_from_start_exact(to_primary_index) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(self.primary_from_start(to_primary_index), cost,);
            } else {
                self.update_primary_from_start(to_primary_index, cost, true);
            }
        }
    }

    pub fn update_additional_secondary_targets(
        &mut self,
        from_secondary_index: AnchorIndex,
        additional_targets: &mut [(SecondaryAnchor, Cost)],
        ts_kind: TsKind,
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_anchor_index, cost) in anchors
            .secondary_anchor_to_index_iter(
                additional_targets.iter().map(|(anchor, _)| *anchor),
                ts_kind,
            )
            .zip(additional_targets.iter().map(|(_, cost)| *cost))
        {
            let Some(to_secondary_index) = to_anchor_index else {
                continue;
            };

            if self.is_secondary_exact(from_secondary_index, to_secondary_index, ts_kind) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(
                    self.secondary(from_secondary_index, to_secondary_index, ts_kind),
                    cost,
                );
            } else {
                self.update_secondary(
                    from_secondary_index,
                    to_secondary_index,
                    ts_kind,
                    cost,
                    true,
                );
            }
        }
    }
}
