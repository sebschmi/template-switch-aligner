use std::fmt::{Display, Formatter, Result, Write};

use crate::costs::cost::Cost;

pub trait IAlignmentType {
    fn is_repeatable(&self) -> bool;

    fn is_repeated(&self, previous: &Self) -> bool;

    fn is_internal(&self) -> bool;
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AlignmentResult<AlignmentType> {
    pub alignment: Vec<(usize, AlignmentType)>,
    pub cost: Cost,
    pub cost_per_base: f64,
    pub duration_seconds: f64,
    pub opened_nodes: usize,
    pub closed_nodes: usize,
    pub suboptimal_opened_nodes: usize,
    pub suboptimal_opened_nodes_ratio: f64,
}

impl<AlignmentType> AlignmentResult<AlignmentType> {
    #[expect(clippy::too_many_arguments)]
    pub fn new(
        alignment: Vec<(usize, AlignmentType)>,
        cost: Cost,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self {
            alignment,
            cost,
            cost_per_base: (cost.as_u64() * 2) as f64 / (reference_length + query_length) as f64,
            duration_seconds,
            opened_nodes,
            closed_nodes,
            suboptimal_opened_nodes,
            suboptimal_opened_nodes_ratio: suboptimal_opened_nodes as f64
                / (opened_nodes - suboptimal_opened_nodes) as f64,
        }
    }
}

impl<AlignmentType: IAlignmentType> AlignmentResult<AlignmentType> {
    pub fn cigar(&self) -> String
    where
        AlignmentType: Display,
    {
        let mut result = String::new();
        self.write_cigar(&mut result).unwrap();
        result
    }

    pub fn write_cigar(&self, writer: &mut impl Write) -> Result
    where
        AlignmentType: Display,
    {
        for (amount, alignment_type) in &self.alignment {
            if alignment_type.is_repeatable() {
                write!(writer, "{amount}{alignment_type}")?;
            } else {
                write!(writer, "{alignment_type}")?;
            }
        }

        Ok(())
    }
}

impl<AlignmentType: Display + IAlignmentType> Display for AlignmentResult<AlignmentType> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        writeln!(f, "Cost: {}", self.cost)?;
        write!(f, "CIGAR: ")?;
        self.write_cigar(f)?;
        writeln!(f)?;
        writeln!(f, "Cost per base: {:.2}", self.cost_per_base)?;
        writeln!(f, "Opened nodes: {}", self.opened_nodes)?;
        writeln!(f, "Closed nodes: {}", self.closed_nodes)?;
        writeln!(
            f,
            "Suboptimal openend nodes: {}",
            self.suboptimal_opened_nodes
        )?;
        writeln!(
            f,
            "Suboptimal openend nodes per optimal opened node: {:.2}",
            self.suboptimal_opened_nodes_ratio
        )?;
        write!(f, "Duration: {:.2}s", self.duration_seconds)?;

        Ok(())
    }
}
