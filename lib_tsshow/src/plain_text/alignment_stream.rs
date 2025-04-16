use std::{collections::VecDeque, iter};

use lib_tsalign::a_star_aligner::{
    alignment_result::alignment::Alignment,
    template_switch_distance::{AlignmentType, TemplateSwitchPrimary},
};

#[derive(Debug, Clone, Default)]
pub struct AlignmentStream {
    stream: VecDeque<(usize, AlignmentType)>,
    length: usize,

    /// The alignment coordinates of the first alignment after the head of the queue.
    head_coordinates: AlignmentCoordinates,
    /// The alignment coordinates of the tailmost element in the queue.
    tail_coordinates: AlignmentCoordinates,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct AlignmentCoordinates {
    reference: usize,
    query: usize,
    template_switch_primary: Option<TemplateSwitchPrimary>,
}

impl AlignmentStream {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn new_with_offset(reference_offset: usize, query_offset: usize) -> Self {
        let mut result = Self::new();
        result.head_coordinates.reference = reference_offset;
        result.head_coordinates.query = query_offset;
        result.tail_coordinates.reference = reference_offset;
        result.tail_coordinates.query = query_offset;
        result
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn stream_iter(&self) -> impl use<'_> + Iterator<Item = (usize, AlignmentType)> {
        self.stream.iter().copied()
    }

    pub fn stream_iter_flat(&self) -> impl use<'_> + Iterator<Item = AlignmentType> {
        self.stream_iter()
            .flat_map(|(multiplicity, alignment_type)| iter::repeat_n(alignment_type, multiplicity))
    }

    fn stream_vec(&self) -> Vec<(usize, AlignmentType)> {
        self.stream_iter().collect()
    }

    pub fn stream_alignment(&self) -> Alignment<AlignmentType> {
        self.stream_vec().into()
    }

    pub fn head_coordinates(&self) -> AlignmentCoordinates {
        self.head_coordinates
    }

    pub fn tail_coordinates(&self) -> AlignmentCoordinates {
        self.tail_coordinates
    }

    pub fn push_until_full(
        &mut self,
        multiplicity: &mut usize,
        alignment_type: AlignmentType,
        requested_length: usize,
    ) {
        let available_length = requested_length - self.length;
        let push_length = *multiplicity * Self::stream_length(alignment_type);

        if available_length >= push_length {
            self.push(*multiplicity, alignment_type);
            *multiplicity = 0;
        } else {
            let push_multiplicity = available_length.div_ceil(Self::stream_length(alignment_type));
            self.push(push_multiplicity, alignment_type);
            *multiplicity -= push_multiplicity;
        }
    }

    pub fn clear(&mut self) {
        self.pop(0);
    }

    pub fn is_full(&self, requested_length: usize) -> bool {
        self.length >= requested_length
    }

    pub fn is_empty(&self) -> bool {
        self.stream_iter().next().is_none()
    }

    pub fn push(&mut self, multiplicity: usize, alignment_type: AlignmentType) {
        self.stream.push_back((multiplicity, alignment_type));
        self.head_coordinates.advance(multiplicity, alignment_type);
        self.length += multiplicity * Self::stream_length(alignment_type);
    }

    /// Pops one unit of length from the tail of the stream.
    pub fn pop_one(&mut self) {
        self.pop(self.length.saturating_sub(1));
    }

    pub fn pop(&mut self, requested_length: usize) {
        while self.length > requested_length {
            let requested_pop_length = self.length - requested_length;
            let (multiplicity, alignment_type) = self.stream.front_mut().unwrap();
            let front_length = *multiplicity * Self::stream_length(*alignment_type);

            if front_length <= requested_pop_length {
                self.tail_coordinates
                    .advance(*multiplicity, *alignment_type);
                self.stream.pop_front();
                self.length -= front_length;
            } else {
                let pop_multiplicity = requested_pop_length / Self::stream_length(*alignment_type);
                self.tail_coordinates
                    .advance(pop_multiplicity, *alignment_type);
                *multiplicity -= pop_multiplicity;
                self.length -= pop_multiplicity * Self::stream_length(*alignment_type);
                break;
            }
        }

        while let Some((multiplicity, alignment_type)) = self.stream.front().copied() {
            if Self::stream_length(alignment_type) == 0 {
                self.tail_coordinates.advance(multiplicity, alignment_type);
                self.stream.pop_front();
            } else {
                break;
            }
        }
    }

    fn stream_length(alignment_type: AlignmentType) -> usize {
        match alignment_type {
            AlignmentType::PrimaryInsertion
            | AlignmentType::PrimaryDeletion
            | AlignmentType::PrimarySubstitution
            | AlignmentType::PrimaryMatch
            | AlignmentType::PrimaryFlankInsertion
            | AlignmentType::PrimaryFlankDeletion
            | AlignmentType::PrimaryFlankSubstitution
            | AlignmentType::PrimaryFlankMatch
            | AlignmentType::SecondaryInsertion
            | AlignmentType::SecondaryDeletion
            | AlignmentType::SecondarySubstitution
            | AlignmentType::SecondaryMatch => 1,
            AlignmentType::TemplateSwitchEntrance { .. }
            | AlignmentType::TemplateSwitchExit { .. }
            | AlignmentType::Root
            | AlignmentType::SecondaryRoot
            | AlignmentType::PrimaryReentry => 0,
            AlignmentType::PrimaryShortcut { .. } => {
                unreachable!("Shortcut alignments are not supported for show")
            }
        }
    }
}

impl AlignmentCoordinates {
    pub fn reference(&self) -> usize {
        self.reference
    }

    pub fn query(&self) -> usize {
        self.query
    }

    fn advance(&mut self, multiplicity: usize, alignment_type: AlignmentType) {
        let (reference_length, query_length) = match alignment_type {
            AlignmentType::PrimaryInsertion | AlignmentType::PrimaryFlankInsertion => (0, 1),
            AlignmentType::PrimaryDeletion | AlignmentType::PrimaryFlankDeletion => (1, 0),
            AlignmentType::PrimarySubstitution
            | AlignmentType::PrimaryMatch
            | AlignmentType::PrimaryFlankSubstitution
            | AlignmentType::PrimaryFlankMatch => (1, 1),
            AlignmentType::TemplateSwitchEntrance { primary, .. } => {
                assert!(
                    self.template_switch_primary.is_none(),
                    "Encountered template switch entrance within template switch"
                );
                self.template_switch_primary = Some(primary);
                (0, 0)
            }
            AlignmentType::SecondaryInsertion
            | AlignmentType::SecondarySubstitution
            | AlignmentType::SecondaryMatch => match self.template_switch_primary.unwrap() {
                TemplateSwitchPrimary::Reference => (1, 0),
                TemplateSwitchPrimary::Query => (0, 1),
            },
            AlignmentType::TemplateSwitchExit { anti_primary_gap } => {
                let Some(template_switch_primary) = self.template_switch_primary.take() else {
                    panic!(
                        "Encountered template switch exit without first encountering a template switch entrance"
                    )
                };
                match template_switch_primary {
                    TemplateSwitchPrimary::Reference => {
                        self.query =
                            usize::try_from(self.query as isize + anti_primary_gap).unwrap()
                    }
                    TemplateSwitchPrimary::Query => {
                        self.reference =
                            usize::try_from(self.reference as isize + anti_primary_gap).unwrap()
                    }
                }
                (0, 0)
            }
            AlignmentType::SecondaryDeletion
            | AlignmentType::Root
            | AlignmentType::SecondaryRoot
            | AlignmentType::PrimaryReentry => (0, 0),
            AlignmentType::PrimaryShortcut { .. } => {
                unreachable!("Shortcut alignments are not supported for show")
            }
        };

        self.reference += multiplicity * reference_length;
        self.query += multiplicity * query_length;
    }
}
