use lib_tsalign::a_star_aligner::template_switch_distance::{
    AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary,
};
use log::debug;

use super::alignment_stream::{AlignmentCoordinates, AlignmentStream};

pub const STREAM_DEFAULT_LENGTH: usize = 20;
pub const STREAM_PADDING: usize = 10;

#[derive(Debug)]
pub struct TSShow<AlignmentType> {
    pub upstream_offset: AlignmentCoordinates,
    pub downstream_limit: AlignmentCoordinates,
    pub sp1_offset: AlignmentCoordinates,
    pub sp2_secondary_offset: usize,
    pub sp3_secondary_offset: usize,
    pub sp4_offset: AlignmentCoordinates,
    pub primary: TemplateSwitchPrimary,
    pub secondary: TemplateSwitchSecondary,
    pub upstream: Vec<(usize, AlignmentType)>,
    pub template_switch: Vec<(usize, AlignmentType)>,
    pub downstream: Vec<(usize, AlignmentType)>,
}

struct AlignmentIterator<'alignment> {
    alignment: &'alignment [(usize, AlignmentType)],
    index: usize,
    current_multiplicity: usize,
}

pub fn parse(alignment: &[(usize, AlignmentType)]) -> Vec<TSShow<AlignmentType>> {
    let mut template_switches = Vec::new();
    let mut stream = AlignmentStream::new();
    let mut alignment = AlignmentIterator::new(alignment);

    while let Some((multiplicity, alignment_type)) = alignment.peek_mut() {
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            template_switches.push(parse_template_switch(&mut alignment, &mut stream));
        } else if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
            panic!("Found template switch exit without matching entrance");
        } else {
            stream.push(*multiplicity, alignment_type);
            alignment.next();
        }
    }

    template_switches
}

fn parse_template_switch(
    alignment: &mut AlignmentIterator,
    stream: &mut AlignmentStream,
) -> TSShow<AlignmentType> {
    let (
        multiplicity,
        alignment_type @ AlignmentType::TemplateSwitchEntrance {
            primary,
            secondary,
            first_offset,
        },
    ) = alignment.peek_mut().unwrap()
    else {
        unreachable!("Function is only called with a template switch entrance")
    };
    debug!("Parsing TS with first_offset {first_offset}");

    let sp1_offset = stream.head_coordinates();
    let upstream = stream.clone();
    let mut template_switch = Vec::new();

    stream.push(*multiplicity, alignment_type);
    alignment.next().unwrap();

    let sp2_secondary_offset: usize = match secondary {
        TemplateSwitchSecondary::Reference => (sp1_offset.reference() as isize + first_offset)
            .try_into()
            .unwrap(),
        TemplateSwitchSecondary::Query => (sp1_offset.query() as isize + first_offset)
            .try_into()
            .unwrap(),
    };
    let mut sp3_secondary_offset = sp2_secondary_offset;

    while let Some((multiplicity, alignment_type)) = alignment.next() {
        debug!("alignment type: {alignment_type}");

        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            panic!("Found template switch entrance within template switch");
        } else if let AlignmentType::TemplateSwitchExit { .. } = alignment_type {
            stream.push(multiplicity, alignment_type);

            let mut upstream = upstream;
            upstream.pop(
                STREAM_DEFAULT_LENGTH.max(
                    sp1_offset
                        .reference()
                        .max(sp1_offset.query())
                        .saturating_sub(sp2_secondary_offset)
                        + STREAM_PADDING,
                ),
            );
            let upstream_offset = upstream.tail_coordinates();
            let upstream = upstream.stream_vec();

            let template_switch = template_switch;
            stream.clear();
            let sp4_offset = stream.head_coordinates();
            let downstream = parse_downstream(
                alignment,
                stream,
                STREAM_DEFAULT_LENGTH.max(sp3_secondary_offset.saturating_sub(
                    sp4_offset.reference().min(sp4_offset.query()) + STREAM_PADDING,
                )),
            );
            let downstream_limit = stream.head_coordinates();

            return TSShow {
                upstream_offset,
                downstream_limit,
                sp1_offset,
                sp2_secondary_offset,
                sp3_secondary_offset,
                sp4_offset,
                primary,
                secondary,
                upstream,
                template_switch,
                downstream,
            };
        } else {
            template_switch.push((multiplicity, alignment_type));
            stream.push(multiplicity, alignment_type);

            if matches!(
                alignment_type,
                AlignmentType::SecondaryDeletion
                    | AlignmentType::SecondarySubstitution
                    | AlignmentType::SecondaryMatch
            ) {
                sp3_secondary_offset -= multiplicity;
                assert!(sp3_secondary_offset < sp2_secondary_offset);
            }
        }
    }

    panic!("Found template switch without exit")
}

fn parse_downstream(
    alignment: &mut AlignmentIterator,
    stream: &mut AlignmentStream,
    requested_length: usize,
) -> Vec<(usize, AlignmentType)> {
    stream.clear();

    while let Some((multiplicity, alignment_type)) = alignment.peek_mut() {
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            break;
        } else if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
            panic!("Found template switch exit without matching entrance");
        } else {
            stream.push_until_full(multiplicity, alignment_type, requested_length);

            if stream.is_full(requested_length) {
                break;
            }
        }
    }

    stream.stream_vec()
}

impl<'alignment> AlignmentIterator<'alignment> {
    fn new(alignment: &'alignment [(usize, AlignmentType)]) -> Self {
        Self {
            alignment,
            index: 0,
            current_multiplicity: alignment
                .first()
                .map(|(multiplicity, _)| {
                    assert!(*multiplicity > 0);
                    *multiplicity
                })
                .unwrap_or(0),
        }
    }

    fn peek_mut(&mut self) -> Option<(&mut usize, AlignmentType)> {
        self.clear_first();

        self.alignment
            .get(self.index)
            .map(|(_, alignment_type)| (&mut self.current_multiplicity, *alignment_type))
    }

    fn next(&mut self) -> Option<(usize, AlignmentType)> {
        self.clear_first();

        self.alignment.get(self.index).map(|(_, alignment_type)| {
            let result = (self.current_multiplicity, *alignment_type);
            self.index += 1;
            self.current_multiplicity = self
                .alignment
                .get(self.index)
                .map(|(multiplicity, _)| {
                    assert!(*multiplicity > 0);
                    *multiplicity
                })
                .unwrap_or(0);
            result
        })
    }

    fn clear_first(&mut self) {
        let max_multiplicity = self
            .alignment
            .get(self.index)
            .map(|(multiplicity, _)| {
                assert!(*multiplicity > 0);
                *multiplicity
            })
            .unwrap_or(0);
        debug_assert!(
            self.current_multiplicity <= max_multiplicity,
            "current_multiplicity > max_multiplicity: {} > {}; index: {}",
            self.current_multiplicity,
            max_multiplicity,
            self.index,
        );

        if self.current_multiplicity == 0 {
            self.index += 1;
            self.current_multiplicity = self
                .alignment
                .get(self.index)
                .map(|(multiplicity, _)| {
                    assert!(*multiplicity > 0);
                    *multiplicity
                })
                .unwrap_or(0);
        }
    }
}
