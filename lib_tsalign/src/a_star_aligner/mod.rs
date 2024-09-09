use std::{collections::HashMap, fmt::Display, hash::Hash, time::Instant};

use alignment_result::AlignmentResult;
use binary_heap_plus::BinaryHeap;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};

use crate::cost::Cost;

pub mod alignment_result;
pub mod gap_affine_edit_distance;
pub mod template_switch_distance;
#[cfg(test)]
mod tests;

/// A node of the alignment graph.
/// The node must implement [`Ord`](std::cmp::Ord), ordering it by its cost, ascending.
/// The graph defined by the node type must be cycle-free.
pub trait AlignmentGraphNode: Sized + Ord {
    /// A unique identifier of the node.
    /// For example, in case of traditional edit distance, this would be the tuple (i, j) indicating which alignment matrix cell this node belongs to.
    type Identifier;

    /// Context required for generating successors.
    type Context;

    /// The type collecting possible alignment types.
    /// For example, in the case of traditional edit distance, this would be a match, a substitution, an insertion or a deletion.
    type AlignmentType;

    /// Create the root node of the alignment graph.
    fn create_root() -> Self;

    /// Generate the successors of this node.
    fn generate_successors<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        context: &Self::Context,
        output: &mut impl Extend<Self>,
    );

    /// Returns the identifier of this node.
    fn identifier(&self) -> &Self::Identifier;

    /// Returns the cost of this node.
    fn cost(&self) -> Cost;

    /// Returns the identifier of the predecessor of this node.
    fn predecessor(&self) -> Option<&Self::Identifier>;

    /// Returns the type of alignment used to generate this node from the predecessor.
    fn predecessor_alignment_type<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        context: &Self::Context,
    ) -> Self::AlignmentType;

    /// Returns true if this node is a target node of the alignment process.
    fn is_target<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        context: &Self::Context,
    ) -> bool;
}

fn a_star_align<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    Node: AlignmentGraphNode,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    context: Node::Context,
) -> AlignmentResult<Node::AlignmentType>
where
    Node::Identifier: Hash + Eq + Clone + Display,
    Node::AlignmentType: Eq,
{
    let start_time = Instant::now();

    let mut closed_list: HashMap<Node::Identifier, Node> = Default::default();
    let mut open_list = BinaryHeap::new_min();
    let mut opened_nodes = 0;
    // Opened nodes that do not have optimal costs.
    let mut suboptimal_opened_nodes = 0;

    open_list.push(Node::create_root());

    let target_node_identifier = loop {
        let Some(node) = open_list.pop() else {
            unreachable!("did not find target node")
        };

        if let Some(previous_visit) = closed_list.get(node.identifier()) {
            // If we have already visited the node, we now must be visiting it with a higher cost.
            debug_assert!(previous_visit.cost() <= node.cost());
            suboptimal_opened_nodes += 1;
            continue;
        }

        if node.is_target(reference, query, &context) {
            let identifier = node.identifier().clone();
            closed_list.insert(node.identifier().clone(), node);
            break identifier;
        }

        let open_nodes_without_new_successors = open_list.len();
        node.generate_successors(reference, query, &context, &mut open_list);
        opened_nodes += open_list.len() - open_nodes_without_new_successors;

        closed_list.insert(node.identifier().clone(), node);
    };

    let mut alignment = Vec::new();
    let mut current_node = closed_list.get(&target_node_identifier).unwrap();
    let cost = current_node.cost();

    loop {
        let alignment_type = current_node.predecessor_alignment_type(reference, query, &context);

        if let Some((count, previous_alignment_type)) = alignment.last_mut() {
            if alignment_type == *previous_alignment_type {
                *count += 1;
            } else {
                alignment.push((1, alignment_type));
            }
        } else {
            alignment.push((1, alignment_type));
        }

        if let Some(predecessor) = current_node.predecessor() {
            current_node = closed_list.get(predecessor).unwrap();
        } else {
            break;
        }
    }

    // Pop root element.
    alignment.pop().unwrap();
    alignment.reverse();

    let end_time = Instant::now();
    let duration = (end_time - start_time).as_secs_f64();

    AlignmentResult::new(
        alignment,
        cost,
        duration,
        opened_nodes,
        closed_list.len(),
        suboptimal_opened_nodes,
        reference.len(),
        query.len(),
    )
}

pub fn gap_affine_edit_distance_a_star_align<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    context: gap_affine_edit_distance::ScoringTable,
) -> AlignmentResult<gap_affine_edit_distance::AlignmentType> {
    a_star_align::<_, _, gap_affine_edit_distance::Node>(reference, query, context)
}

pub fn template_switch_distance_a_star_align<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    context: template_switch_distance::ScoringTable,
) -> AlignmentResult<template_switch_distance::AlignmentType> {
    a_star_align::<_, _, template_switch_distance::Node>(reference, query, context)
}
