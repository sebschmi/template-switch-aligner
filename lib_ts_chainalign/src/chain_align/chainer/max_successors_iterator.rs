use generic_a_star::cost::AStarCost;

use crate::anchors::index::AnchorIndex;

pub struct MaxSuccessorsIterator<Iter: Iterator<Item = (AnchorIndex, Cost)>, Cost> {
    iter: Iter,
    count: usize,
    max_count: usize,
    initial_cost: Cost,
    next_cost: Option<Cost>,
}

impl<Iter: Iterator<Item = (AnchorIndex, Cost)>, Cost: AStarCost>
    MaxSuccessorsIterator<Iter, Cost>
{
    pub fn new(iter: Iter, max_count: usize) -> Self {
        Self {
            iter,
            count: 0,
            max_count,
            initial_cost: Cost::zero(),
            next_cost: None,
        }
    }

    pub fn successor_count(&self) -> usize {
        self.count
    }

    pub fn next_cost(&self) -> Option<Cost> {
        self.next_cost
    }
}

impl<Iter: Iterator<Item = (AnchorIndex, Cost)>, Cost: AStarCost> Iterator
    for MaxSuccessorsIterator<Iter, Cost>
{
    type Item = (AnchorIndex, Cost);

    fn next(&mut self) -> Option<Self::Item> {
        let (anchor_index, cost) = self.iter.next()?;

        if cost == Cost::max_value() {
            return None;
        }

        if self.count == 0 {
            self.count += 1;
            self.initial_cost = cost;
            Some((anchor_index, cost))
        } else if self.count < self.max_count || cost == self.initial_cost {
            self.count += 1;
            Some((anchor_index, cost))
        } else {
            self.next_cost = Some(cost);
            None
        }
    }
}
