use std::collections::VecDeque;

use crate::{a_star_aligner::template_switch_distance::Context, cost::Cost};

use super::AlignmentStrategy;

pub trait FlankCostStrategy: AlignmentStrategy {
    fn current_left_flank_cost(&self) -> Cost;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoFlankCostStrategy;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SequenceIdentityFlankCostStrategy {
    // TODO have left and right buffers
    buffer: VecDeque<Cost>,
    buffer_sum: Cost,
}

impl FlankCostStrategy for NoFlankCostStrategy {
    fn current_left_flank_cost(&self) -> Cost {
        Cost::ZERO
    }
}

impl FlankCostStrategy for SequenceIdentityFlankCostStrategy {
    fn current_left_flank_cost(&self) -> Cost {
        self.buffer_sum
    }
}

impl AlignmentStrategy for NoFlankCostStrategy {
    fn create_root(_context: &Context) -> Self {
        Self
    }

    fn generate_successor(&self, _context: &Context) -> Self {
        *self
    }
}

impl AlignmentStrategy for SequenceIdentityFlankCostStrategy {
    fn create_root(context: &Context) -> Self {
        assert!(context.flank_length > 0, "Template switch scoring flank length cannot be zero for the chosen flank cost strategy");

        Self {
            buffer: VecDeque::with_capacity(context.flank_length),
            buffer_sum: Cost::ZERO,
        }
    }

    fn generate_successor(&self, context: &Context) -> Self {
        let mut successor = self.clone();
        if successor.buffer.len() == context.flank_length {
            successor.buffer_sum -= successor.buffer.pop_front().unwrap();
        }

        todo!("add new cost to back of buffer")
    }
}
