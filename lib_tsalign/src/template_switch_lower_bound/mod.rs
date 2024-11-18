use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::Cost;
use ndarray::Array2;

use crate::config::TemplateSwitchConfig;

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchLowerBoundMatrix {
    matrix: Array2<Cost>,
    min_distance_between_two_template_switches: usize,
    template_switch_min_length: isize,
    template_switch_max_length: isize,
}

impl TemplateSwitchLowerBoundMatrix {
    pub fn new<AlphabetType: Alphabet>(config: &TemplateSwitchConfig<AlphabetType>) -> Self {
        let min_distance_between_two_template_switches =
            usize::try_from(config.left_flank_length + config.right_flank_length).unwrap();
        let template_switch_min_length =
            isize::try_from(config.length_costs.minimum_finite_input().unwrap()).unwrap()
                + config
                    .length_difference_costs
                    .minimum_finite_input()
                    .unwrap()
                    .min(0);
        let template_switch_max_length =
            isize::try_from(config.length_costs.maximum_finite_input().unwrap()).unwrap()
                + config
                    .length_difference_costs
                    .maximum_finite_input()
                    .unwrap();
        let template_switch_length_range =
            usize::try_from(template_switch_min_length.abs() + template_switch_max_length.abs())
                .unwrap();
        let mut matrix = Array2::zeros((
            template_switch_length_range + 1,
            template_switch_length_range + 1,
        ));

        for (x, y) in (0..template_switch_length_range + 1)
            .flat_map(|x| (x..template_switch_length_range + 1).map(move |y| (x, y)))
        {
            let dx = isize::try_from(x).unwrap() + template_switch_min_length;
            let dy = isize::try_from(y).unwrap() + template_switch_min_length;

            todo!()
        }

        Self {
            matrix,
            min_distance_between_two_template_switches,
            template_switch_min_length,
            template_switch_max_length,
        }
    }

    pub fn min_distance_between_two_template_switches(&self) -> usize {
        self.min_distance_between_two_template_switches
    }

    pub fn template_switch_min_length(&self) -> isize {
        self.template_switch_min_length
    }

    pub fn template_switch_max_length(&self) -> isize {
        self.template_switch_max_length
    }
}
