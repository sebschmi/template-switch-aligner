#![forbid(clippy::mod_module_files)]

use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    hash::Hash,
};

use binary_heap_plus::{BinaryHeap, MinComparator};
use cost::Cost;
use deterministic_default_hasher::DeterministicDefaultHasher;
use reset::Reset;

pub mod cost;
pub mod reset;

/// A node of the A* graph.
/// The node must implement [`Ord`](std::cmp::Ord), ordering it by its cost, ascending.
/// The graph defined by the node type must be cycle-free.
pub trait AStarNode: Sized + Ord + Debug + Display {
    /// A unique identifier of the node.
    ///
    /// For example, in case of traditional edit distance, this would be the tuple (i, j) indicating which alignment matrix cell this node belongs to.
    type Identifier: Debug + Clone + Eq + Hash;

    /// The type collecting possible edge types.
    ///
    /// These are used when backtracking a solution.
    type EdgeType: Debug;

    /// Returns the identifier of this node.
    fn identifier(&self) -> &Self::Identifier;

    /// Returns the cost of this node.
    ///
    /// This is the cost measured from the root node, and does NOT include the A* lower bound.
    fn cost(&self) -> Cost;

    /// Returns the A* lower bound of this node.
    fn a_star_lower_bound(&self) -> Cost;

    /// Returns the identifier of the predecessor of this node.
    fn predecessor(&self) -> Option<&Self::Identifier>;

    /// Returns the edge type used to reach this node from the predecessor, or `None` if this is a root node.
    fn predecessor_edge_type(&self) -> Option<Self::EdgeType>;
}

pub trait AStarContext: Reset {
    /// The node type used by the A* algorithm.
    type Node: AStarNode;

    /// Create the root node of the A* graph.
    fn create_root(&self) -> Self::Node;

    /// Generate the successors of this node.
    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>);

    /// Returns true if this node is a target node of the A* graph.
    fn is_target(&self, node: &Self::Node) -> bool;

    /// Returns the maximum cost that the target node is allowed to have.
    /// If no target is found with this cost or lower, then [`AStarResult::NoTarget`] is returned.
    fn max_cost(&self) -> Option<Cost>;
}

#[derive(Debug, Default)]
pub struct AStarPerformanceCounters {
    pub opened_nodes: usize,
    /// Opened nodes that do not have optimal costs.
    pub suboptimal_opened_nodes: usize,
    pub closed_nodes: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub enum AStarState<NodeIdentifier> {
    /// The algorithm was just created or reset.
    Empty,
    /// The algorithm was just initialised.
    Init,
    /// The algorithm is searching for a target node.
    Searching,
    /// The algorithm terminated.
    Terminated { result: AStarResult<NodeIdentifier> },
}

#[derive(Debug)]
pub struct AStar<Context: AStarContext> {
    state: AStarState<<Context::Node as AStarNode>::Identifier>,
    context: Context,
    closed_list: HashMap<
        <Context::Node as AStarNode>::Identifier,
        Context::Node,
        DeterministicDefaultHasher,
    >,
    open_list: BinaryHeap<Context::Node, MinComparator>,
    performance_counters: AStarPerformanceCounters,
}

#[derive(Debug)]
pub struct AStarBuffers<NodeIdentifier, Node> {
    closed_list: HashMap<NodeIdentifier, Node, DeterministicDefaultHasher>,
    open_list: BinaryHeap<Node, MinComparator>,
}

#[derive(Debug, PartialEq, Eq)]
pub enum AStarResult<NodeIdentifier> {
    /// The algorithm has found a target node.
    FoundTarget {
        identifier: NodeIdentifier,
        cost: Cost,
    },
    /// The algorithm terminated, but did not find a target.
    NoTarget { max_cost: Option<Cost> },
}

struct BacktrackingIterator<'a_star, Context: AStarContext> {
    a_star: &'a_star AStar<Context>,
    current: <Context::Node as AStarNode>::Identifier,
}

struct BacktrackingIteratorWithCost<'a_star, Context: AStarContext> {
    a_star: &'a_star AStar<Context>,
    current: <Context::Node as AStarNode>::Identifier,
}

impl<Context: AStarContext> AStar<Context> {
    pub fn new(context: Context) -> Self {
        Self {
            state: AStarState::Empty,
            context,
            closed_list: Default::default(),
            open_list: BinaryHeap::new_min(),
            performance_counters: Default::default(),
        }
    }

    pub fn new_with_buffers(
        context: Context,
        mut buffers: AStarBuffers<<Context::Node as AStarNode>::Identifier, Context::Node>,
    ) -> Self {
        buffers.closed_list.clear();
        buffers.open_list.clear();
        Self {
            state: AStarState::Empty,
            context,
            closed_list: buffers.closed_list,
            open_list: buffers.open_list,
            performance_counters: Default::default(),
        }
    }

    pub fn state(&self) -> &AStarState<<Context::Node as AStarNode>::Identifier> {
        &self.state
    }

    pub fn context(&self) -> &Context {
        &self.context
    }

    pub fn into_context(self) -> Context {
        self.context
    }

    pub fn into_buffers(
        self,
    ) -> AStarBuffers<<Context::Node as AStarNode>::Identifier, Context::Node> {
        AStarBuffers {
            closed_list: self.closed_list,
            open_list: self.open_list,
        }
    }

    pub fn closed_node(
        &self,
        node_identifier: &<Context::Node as AStarNode>::Identifier,
    ) -> Option<&Context::Node> {
        self.closed_list.get(node_identifier)
    }

    pub fn performance_counters(&self) -> &AStarPerformanceCounters {
        &self.performance_counters
    }

    pub fn reset(&mut self) {
        self.state = AStarState::Empty;
        self.context.reset();
        self.closed_list.clear();
        self.open_list.clear();
        self.performance_counters = Default::default();
    }

    pub fn initialise(&mut self) {
        self.initialise_with(|context| context.create_root());
    }

    pub fn initialise_with(&mut self, node: impl FnOnce(&Context) -> Context::Node) {
        assert_eq!(self.state, AStarState::Empty);

        self.state = AStarState::Init;
        self.open_list.push(node(&self.context));
    }

    pub fn search(&mut self) -> AStarResult<<Context::Node as AStarNode>::Identifier> {
        self.search_until(|context, node| context.is_target(node))
    }

    pub fn search_until(
        &mut self,
        mut is_target: impl FnMut(&Context, &Context::Node) -> bool,
    ) -> AStarResult<<Context::Node as AStarNode>::Identifier> {
        assert!(matches!(
            self.state,
            AStarState::Init | AStarState::Searching | AStarState::Terminated { .. }
        ));

        let max_cost = self.context.max_cost();
        if self.open_list.is_empty() {
            return AStarResult::NoTarget { max_cost };
        }

        self.state = AStarState::Searching;

        let mut last_node = None;

        let target_identifier = loop {
            let Some(node) = self.open_list.pop() else {
                if last_node.is_none() {
                    unreachable!("Open list was empty.");
                };
                self.state = AStarState::Terminated {
                    result: AStarResult::NoTarget { max_cost },
                };
                return AStarResult::NoTarget { max_cost };
            };

            if node.cost() > max_cost.unwrap_or(Cost::MAX) {
                self.state = AStarState::Terminated {
                    result: AStarResult::NoTarget { max_cost },
                };
                return AStarResult::NoTarget { max_cost };
            }

            last_node = Some(node.identifier().clone());

            if let Some(previous_visit) = self.closed_list.get(node.identifier()) {
                // If we have already visited the node, we now must be visiting it with a higher cost.
                debug_assert!(
                    previous_visit.cost() + previous_visit.a_star_lower_bound()
                        <= node.cost() + node.a_star_lower_bound(),
                    "{}",
                    {
                        use std::fmt::Write;
                        let mut previous_visit = previous_visit;
                        let mut node = &node;
                        let mut out = String::new();

                        writeln!(out, "previous_visit:").unwrap();
                        while let Some(predecessor) = previous_visit.predecessor() {
                            writeln!(out, "{previous_visit}").unwrap();
                            previous_visit = self.closed_list.get(predecessor).unwrap();
                        }

                        writeln!(out, "\nnode:").unwrap();
                        while let Some(predecessor) = node.predecessor() {
                            writeln!(out, "{node}").unwrap();
                            node = self.closed_list.get(predecessor).unwrap();
                        }

                        out
                    }
                );
                self.performance_counters.suboptimal_opened_nodes += 1;
                continue;
            }

            let open_nodes_without_new_successors = self.open_list.len();
            self.context.generate_successors(&node, &mut self.open_list);
            self.performance_counters.opened_nodes +=
                self.open_list.len() - open_nodes_without_new_successors;

            if is_target(&self.context, &node) {
                let identifier = node.identifier().clone();
                let previous_visit = self.closed_list.insert(node.identifier().clone(), node);
                self.performance_counters.closed_nodes += 1;
                debug_assert!(previous_visit.is_none());
                break identifier;
            }

            let previous_visit = self.closed_list.insert(node.identifier().clone(), node);
            self.performance_counters.closed_nodes += 1;
            debug_assert!(previous_visit.is_none());
        };

        let cost = self.closed_list.get(&target_identifier).unwrap().cost();
        self.state = AStarState::Terminated {
            result: AStarResult::FoundTarget {
                identifier: target_identifier.clone(),
                cost,
            },
        };
        AStarResult::FoundTarget {
            identifier: target_identifier,
            cost,
        }
    }

    pub fn backtrack(
        &self,
    ) -> impl use<'_, Context> + Iterator<Item = <Context::Node as AStarNode>::EdgeType> {
        let AStarState::Terminated {
            result: AStarResult::FoundTarget { identifier, .. },
        } = &self.state
        else {
            panic!("Cannot backtrack since no target was found.")
        };

        self.backtrack_from(identifier).unwrap()
    }

    /// Backtrack from the target node to a root node.
    ///
    /// The elements of the iterator are a pair of an edge and the cost of the node that is reached by the edge.
    /// The cost of the first node is never returned.
    pub fn backtrack_with_costs(
        &self,
    ) -> impl use<'_, Context> + Iterator<Item = (<Context::Node as AStarNode>::EdgeType, Cost)>
    {
        let AStarState::Terminated {
            result: AStarResult::FoundTarget { identifier, .. },
        } = &self.state
        else {
            panic!("Cannot backtrack since no target was found.")
        };

        self.backtrack_with_costs_from(identifier).unwrap()
    }

    pub fn backtrack_from(
        &self,
        identifier: &<Context::Node as AStarNode>::Identifier,
    ) -> Option<impl use<'_, Context> + Iterator<Item = <Context::Node as AStarNode>::EdgeType>>
    {
        if self.closed_list.contains_key(identifier) {
            Some(BacktrackingIterator {
                a_star: self,
                current: identifier.clone(),
            })
        } else {
            None
        }
    }

    /// Backtrack from a node to a root node.
    ///
    /// The elements of the iterator are a pair of an edge and the cost of the node that is reached by the edge.
    /// The cost of the first node is never returned.
    pub fn backtrack_with_costs_from(
        &self,
        identifier: &<Context::Node as AStarNode>::Identifier,
    ) -> Option<
        impl use<'_, Context> + Iterator<Item = (<Context::Node as AStarNode>::EdgeType, Cost)>,
    > {
        if self.closed_list.contains_key(identifier) {
            Some(BacktrackingIteratorWithCost {
                a_star: self,
                current: identifier.clone(),
            })
        } else {
            None
        }
    }
}

impl<Context: AStarContext> Iterator for BacktrackingIterator<'_, Context> {
    type Item = <Context::Node as AStarNode>::EdgeType;

    fn next(&mut self) -> Option<Self::Item> {
        let current = self.a_star.closed_list.get(&self.current).unwrap();

        if let Some(predecessor) = current.predecessor().cloned() {
            let predecessor_edge_type = current.predecessor_edge_type().unwrap();
            self.current = predecessor;
            Some(predecessor_edge_type)
        } else {
            None
        }
    }
}

impl<Context: AStarContext> Iterator for BacktrackingIteratorWithCost<'_, Context> {
    type Item = (<Context::Node as AStarNode>::EdgeType, Cost);

    fn next(&mut self) -> Option<Self::Item> {
        let current = self.a_star.closed_list.get(&self.current).unwrap();
        let cost = current.cost();

        if let Some(predecessor) = current.predecessor().cloned() {
            let predecessor_edge_type = current.predecessor_edge_type().unwrap();
            self.current = predecessor;
            Some((predecessor_edge_type, cost))
        } else {
            None
        }
    }
}

impl<NodeIdentifier, Node: Ord> Default for AStarBuffers<NodeIdentifier, Node> {
    fn default() -> Self {
        Self {
            closed_list: Default::default(),
            open_list: BinaryHeap::new_min(),
        }
    }
}
