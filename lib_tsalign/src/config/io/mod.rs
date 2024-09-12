use compact_genome::interface::alphabet::Alphabet;

use crate::{costs::gap_affine::GapAffineAlignmentCostTable, error::Error};

use super::TemplateSwitchConfig;

impl<AlphabetType: Alphabet> TemplateSwitchConfig<AlphabetType> {
    pub fn read_plain(reader: impl std::io::Read) -> crate::error::Result<Self> {
        let cost_tables = GapAffineAlignmentCostTable::read_plain_multi(reader)?;

        let mut keys: Vec<_> = cost_tables.keys().cloned().collect();
        let mut expected_keys: Vec<_> = [
            "Primary Edit Costs",
            "Secondary Edit Costs",
            "Left Flank Edit Costs",
            "Right Flank Edit Costs",
        ]
        .into_iter()
        .map(ToString::to_string)
        .collect();

        keys.sort_unstable();
        expected_keys.sort_unstable();

        if keys != expected_keys {
            return Err(Error::WrongCostTableNames {
                actual: keys,
                expected: expected_keys,
            });
        }

        #[expect(unreachable_code)]
        Ok(Self {
            left_flank_length: todo!(),
            right_flank_length: todo!(),

            primary_edit_costs: cost_tables.remove("Primary Edit Costs").unwrap(),
            secondary_edit_costs: cost_tables.remove("Secondary Edit Costs").unwrap(),
            left_flank_edit_costs: cost_tables.remove("Left Flank Edit Costs").unwrap(),
            right_flank_edit_costs: cost_tables.remove("Right Flank Edit Costs").unwrap(),

            offset1_costs: todo!(),
        })
    }
}
