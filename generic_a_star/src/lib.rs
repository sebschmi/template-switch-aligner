#![forbid(clippy::mod_module_files)]

use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    hash::Hash,
};

use binary_heap_plus::{BinaryHeap, MinComparator};
use cost::AStarCost;
use deterministic_default_hasher::DeterministicDefaultHasher;
use extend_map::ExtendFilter;
use num_traits::Bounded;
use reset::Reset;

pub mod cost;
pub mod reset;

/// A node of the A* graph.
/// The node must implement [`Ord`], ordering it by its cost plus A* cost, ascending.
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

    type Cost: AStarCost;

    /// Returns the identifier of this node.
    fn identifier(&self) -> &Self::Identifier;

    /// Returns the cost of this node.
    ///
    /// This is the cost measured from the root node, and does NOT include the A* lower bound.
    fn cost(&self) -> Self::Cost;

    /// Returns the A* lower bound of this node.
    fn a_star_lower_bound(&self) -> Self::Cost;

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
    ///
    /// If no target is found with this cost or lower, then [`AStarResult::ExceededCostLimit`] is returned.
    fn cost_limit(&self) -> Option<<Self::Node as AStarNode>::Cost>;

    /// An approximate memory limit for the aligner in bytes.
    ///
    /// If it is exceeded, then [`AStarResult::ExceededMemoryLimit`] is returned
    fn memory_limit(&self) -> Option<usize>;

    /// Returns true if the nodes are generated in a label-setting manner.
    ///
    /// Label setting means that once a node has been closed, it will never be opened at a smaller cost.
    /// On the contrary, if this is set to false, then nodes are allowed to be closed multiple times.
    /// This results in a worse performance, but allows for example to handle negative costs.
    ///
    /// This method returns `true` in its default implementation.
    fn is_label_setting(&self) -> bool {
        true
    }
}

#[derive(Debug, Default)]
pub struct AStarPerformanceCounters {
    pub opened_nodes: usize,
    /// Opened nodes that do not have optimal costs.
    pub suboptimal_opened_nodes: usize,
    pub closed_nodes: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub enum AStarState<NodeIdentifier, Cost> {
    /// The algorithm was just created or reset.
    Empty,
    /// The algorithm was just initialised.
    Init,
    /// The algorithm is searching for a target node.
    Searching,
    /// The algorithm terminated.
    Terminated {
        result: AStarResult<NodeIdentifier, Cost>,
    },
}

#[derive(Debug)]
pub struct AStar<Context: AStarContext> {
    state: AStarState<<Context::Node as AStarNode>::Identifier, <Context::Node as AStarNode>::Cost>,
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

#[derive(Debug, Clone, Ord, PartialOrd, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "astar_result_type"))]
pub enum AStarResult<NodeIdentifier, Cost> {
    /// The algorithm has found a target node.
    FoundTarget {
        #[cfg_attr(feature = "serde", serde(skip))]
        identifier: NodeIdentifier,
        cost: Cost,
    },

    /// The algorithm terminated before finding a target because the cost limit was reached.
    ExceededCostLimit { cost_limit: Cost },

    /// The algorithm termianted before finding a target because the memory limit was reached.
    ExceededMemoryLimit {
        /// The maximum cost reached before reaching the memory limit.
        max_cost: Cost,
    },

    /// The algorithm terminated, but did not find a target.
    NoTarget,
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

    pub fn state(
        &self,
    ) -> &AStarState<<Context::Node as AStarNode>::Identifier, <Context::Node as AStarNode>::Cost>
    {
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

    pub fn search(
        &mut self,
    ) -> AStarResult<<Context::Node as AStarNode>::Identifier, <Context::Node as AStarNode>::Cost>
    {
        self.search_until(|context, node| context.is_target(node))
    }

    pub fn search_until(
        &mut self,
        mut is_target: impl FnMut(&Context, &Context::Node) -> bool,
    ) -> AStarResult<<Context::Node as AStarNode>::Identifier, <Context::Node as AStarNode>::Cost>
    {
        assert!(matches!(
            self.state,
            AStarState::Init | AStarState::Searching | AStarState::Terminated { .. }
        ));

        let cost_limit = self
            .context
            .cost_limit()
            .unwrap_or(<Context::Node as AStarNode>::Cost::max_value());
        let mut applied_cost_limit = false;
        let memory_limit = self.context.memory_limit().unwrap_or(usize::MAX);
        // The factor of 2.3 is determined empirically.
        let node_count_limit =
            (memory_limit as f64 / std::mem::size_of::<Context::Node>() as f64 / 2.3).round()
                as usize;

        if self.open_list.is_empty() {
            return AStarResult::NoTarget;
        }

        self.state = AStarState::Searching;

        let mut last_node = None;

        let target_identifier = loop {
            let Some(node) = self.open_list.pop() else {
                if last_node.is_none() {
                    unreachable!("Open list was empty.");
                };
                if applied_cost_limit {
                    self.state = AStarState::Terminated {
                        result: AStarResult::ExceededCostLimit { cost_limit },
                    };
                    return AStarResult::ExceededCostLimit { cost_limit };
                } else {
                    self.state = AStarState::Terminated {
                        result: AStarResult::NoTarget,
                    };
                    return AStarResult::NoTarget;
                }
            };

            // Nodes are ordered by cost plus lower bound.
            if node.cost() + node.a_star_lower_bound() > cost_limit {
                self.state = AStarState::Terminated {
                    result: AStarResult::ExceededCostLimit { cost_limit },
                };
                return AStarResult::ExceededCostLimit { cost_limit };
            }

            if self.closed_list.len() + self.open_list.len() > node_count_limit {
                self.state = AStarState::Terminated {
                    result: AStarResult::ExceededMemoryLimit {
                        max_cost: node.cost(),
                    },
                };
                return AStarResult::ExceededMemoryLimit {
                    max_cost: node.cost(),
                };
            }

            last_node = Some(node.identifier().clone());

            if let Some(previous_visit) = self.closed_list.get(node.identifier()) {
                if self.context.is_label_setting() {
                    // In label-setting mode, if we have already visited the node, we now must be visiting it with a higher cost.
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
                }

                self.performance_counters.suboptimal_opened_nodes += 1;
                continue;
            }

            let open_nodes_without_new_successors = self.open_list.len();
            self.context.generate_successors(
                &node,
                &mut ExtendFilter::new(&mut self.open_list, |node| {
                    let result = node.cost() + node.a_star_lower_bound() <= cost_limit;
                    applied_cost_limit = applied_cost_limit || !result;
                    result
                }),
            );
            self.performance_counters.opened_nodes +=
                self.open_list.len() - open_nodes_without_new_successors;

            if is_target(&self.context, &node) {
                let identifier = node.identifier().clone();
                let previous_visit = self.closed_list.insert(node.identifier().clone(), node);
                self.performance_counters.closed_nodes += 1;
                debug_assert!(previous_visit.is_none() || !self.context.is_label_setting());
                break identifier;
            }

            let previous_visit = self.closed_list.insert(node.identifier().clone(), node);
            self.performance_counters.closed_nodes += 1;
            debug_assert!(previous_visit.is_none() || !self.context.is_label_setting());
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
    ) -> impl use<'_, Context>
    + Iterator<
        Item = (
            <Context::Node as AStarNode>::EdgeType,
            <Context::Node as AStarNode>::Cost,
        ),
    > {
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
    #[allow(clippy::type_complexity)]
    pub fn backtrack_with_costs_from(
        &self,
        identifier: &<Context::Node as AStarNode>::Identifier,
    ) -> Option<
        impl use<'_, Context>
        + Iterator<
            Item = (
                <Context::Node as AStarNode>::EdgeType,
                <Context::Node as AStarNode>::Cost,
            ),
        >,
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

impl<NodeIdentifier, Cost: Copy> AStarResult<NodeIdentifier, Cost> {
    /// Returns the maximum cost of closed nodes reached during alignment.
    ///
    /// **Panics** if `self` is [`AStarResult::NoTarget`].
    pub fn cost(&self) -> Cost {
        match self {
            Self::FoundTarget { cost, .. } => *cost,
            Self::ExceededCostLimit { cost_limit } => *cost_limit,
            Self::ExceededMemoryLimit { max_cost } => *max_cost,
            Self::NoTarget => panic!("AStarResult has no costs"),
        }
    }

    pub fn without_node_identifier(&self) -> AStarResult<(), Cost> {
        match *self {
            Self::FoundTarget { cost, .. } => AStarResult::FoundTarget {
                identifier: (),
                cost,
            },
            Self::ExceededCostLimit { cost_limit } => AStarResult::ExceededCostLimit { cost_limit },
            Self::ExceededMemoryLimit { max_cost } => AStarResult::ExceededMemoryLimit { max_cost },
            Self::NoTarget => AStarResult::NoTarget,
        }
    }
}

impl<NodeIdentifier: Clone, Cost> AStarResult<NodeIdentifier, Cost> {
    pub fn transform_cost<TargetCost>(
        &self,
        transform: impl Fn(&Cost) -> TargetCost,
    ) -> AStarResult<NodeIdentifier, TargetCost> {
        match self {
            AStarResult::FoundTarget { identifier, cost } => AStarResult::FoundTarget {
                identifier: identifier.clone(),
                cost: transform(cost),
            },
            AStarResult::ExceededCostLimit { cost_limit } => AStarResult::ExceededCostLimit {
                cost_limit: transform(cost_limit),
            },
            AStarResult::ExceededMemoryLimit { max_cost } => AStarResult::ExceededMemoryLimit {
                max_cost: transform(max_cost),
            },
            AStarResult::NoTarget => AStarResult::NoTarget,
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
    type Item = (
        <Context::Node as AStarNode>::EdgeType,
        <Context::Node as AStarNode>::Cost,
    );

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

impl<NodeIdentifier, Cost: Display> Display for AStarResult<NodeIdentifier, Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AStarResult::FoundTarget { cost, .. } => write!(f, "Reached target with cost {cost}"),
            AStarResult::ExceededCostLimit { cost_limit } => {
                write!(f, "Exceeded cost limit of {cost_limit}")
            }
            AStarResult::ExceededMemoryLimit { max_cost } => write!(
                f,
                "Exceeded memory limit, but reached a maximum cost of {max_cost}"
            ),
            AStarResult::NoTarget => write!(f, "Found no target"),
        }
    }
}

impl<NodeIdentifier, Cost> Default for AStarResult<NodeIdentifier, Cost> {
    fn default() -> Self {
        Self::NoTarget
    }
}
