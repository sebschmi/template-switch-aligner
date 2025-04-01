use std::collections::{BTreeMap, HashMap};

use builder::MultipairAlignmentImplementation;
use coordinates::{ArrangementColumn, ArrangementRow, SourceColumn};
use tagged_vec::TaggedVec;

mod alignment;
mod builder;
pub mod coordinates;
mod sequence;

pub struct MultipairAlignmentArrangement<RowKey = String> {
    alignment: MultipairAlignment<RowKey>,
    arrangement: HashMap<ArrangementRow, TaggedVec<SourceColumn, ArrangementColumn>>,
}

pub struct MultipairAlignment<RowKey = String> {
    implementation: MultipairAlignmentImplementation,
    row_keys_to_arrangement_rows: BTreeMap<RowKey, ArrangementRow>,
}
