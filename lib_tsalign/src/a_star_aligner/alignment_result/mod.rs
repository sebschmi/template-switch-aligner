use std::fmt::{Display, Formatter, Result, Write};

use crate::cost::Cost;

pub struct AlignmentResult<AlignmentType> {
    pub alignment: Vec<(usize, AlignmentType)>,
    pub cost: Cost,
    pub cost_per_base: f64,
    pub opened_nodes: usize,
    pub closed_nodes: usize,
}

impl<AlignmentType> AlignmentResult<AlignmentType> {
    pub fn new(
        alignment: Vec<(usize, AlignmentType)>,
        cost: Cost,
        opened_nodes: usize,
        closed_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self {
            alignment,
            cost,
            cost_per_base: (cost.as_u64() * 2) as f64 / (reference_length + query_length) as f64,
            opened_nodes,
            closed_nodes,
        }
    }

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
            write!(writer, "{amount}{alignment_type}")?;
        }

        Ok(())
    }
}

impl<AlignmentType: Display> Display for AlignmentResult<AlignmentType> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        writeln!(f, "Cost: {}", self.cost)?;
        writeln!(f, "Cost per base: {:.2}", self.cost_per_base)?;
        writeln!(f, "Opened nodes: {}", self.opened_nodes)?;
        writeln!(f, "Closed nodes: {}", self.closed_nodes)?;
        write!(f, "CIGAR: ")?;
        self.write_cigar(f)?;

        Ok(())
    }
}
