use lib_tsalign::a_star_aligner::template_switch_distance::{
    AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary,
};

use super::alignment_stream::{AlignmentCoordinates, AlignmentStream};

#[derive(Debug)]
pub struct TSShow<AlignmentType> {
    pub upstream_offset: AlignmentCoordinates,
    pub downstream_limit: AlignmentCoordinates,
    pub sp1_offset: AlignmentCoordinates,
    pub sp2_offset: usize,
    pub sp4_offset: AlignmentCoordinates,
    pub primary: TemplateSwitchPrimary,
    pub secondary: TemplateSwitchSecondary,
    pub first_offset: isize,
    pub length_difference: isize,
    pub upstream: Vec<(usize, AlignmentType)>,
    pub template_switch: Vec<(usize, AlignmentType)>,
    pub downstream: Vec<(usize, AlignmentType)>,
}

pub fn parse(mut alignment: &[(usize, AlignmentType)]) -> Vec<TSShow<AlignmentType>> {
    let mut template_switches = Vec::new();
    let mut stream = AlignmentStream::new(20);

    while let Some((multiplicity, alignment_type)) = alignment.first().copied() {
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            template_switches.push(parse_template_switch(&mut alignment, &mut stream));
        } else if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
            panic!("Found template switch exit without matching entrance");
        } else {
            stream.push(multiplicity, alignment_type);
            alignment = &alignment[1..];
        }
    }

    template_switches
}

fn parse_template_switch(
    alignment: &mut &[(usize, AlignmentType)],
    stream: &mut AlignmentStream,
) -> TSShow<AlignmentType> {
    let upstream_offset = stream.tail_coordinates();
    let sp1_offset = stream.head_coordinates();
    let upstream = stream.stream_vec();
    let mut template_switch = Vec::new();

    let (
        multiplicity,
        alignment_type @ AlignmentType::TemplateSwitchEntrance {
            primary,
            secondary,
            first_offset,
        },
    ) = alignment.first().copied().unwrap()
    else {
        unreachable!("Function is only called with a template switch entrance")
    };
    template_switch.push((multiplicity, alignment_type));
    stream.push(multiplicity, alignment_type);
    *alignment = &alignment[1..];

    let sp2_offset = match secondary {
        TemplateSwitchSecondary::Reference => (sp1_offset.reference() as isize + first_offset)
            .try_into()
            .unwrap(),
        TemplateSwitchSecondary::Query => (sp1_offset.query() as isize + first_offset)
            .try_into()
            .unwrap(),
    };

    while let Some((multiplicity, alignment_type)) = alignment.first().copied() {
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            panic!("Found template switch entrance within template switch");
        } else if let AlignmentType::TemplateSwitchExit { length_difference } = alignment_type {
            template_switch.push((multiplicity, alignment_type));
            stream.push(multiplicity, alignment_type);
            *alignment = &alignment[1..];

            let template_switch = template_switch;
            stream.clear();
            let sp4_offset = stream.head_coordinates();
            let downstream = parse_downstream(alignment, stream);
            let downstream_limit = stream.head_coordinates();

            return TSShow {
                upstream_offset,
                downstream_limit,
                sp1_offset,
                sp2_offset,
                sp4_offset,
                primary,
                secondary,
                first_offset,
                length_difference,
                upstream,
                template_switch,
                downstream,
            };
        } else {
            template_switch.push((multiplicity, alignment_type));
            stream.push(multiplicity, alignment_type);
            *alignment = &alignment[1..];
        }
    }

    panic!("Found template switch without exit")
}

fn parse_downstream(
    alignment: &mut &[(usize, AlignmentType)],
    stream: &mut AlignmentStream,
) -> Vec<(usize, AlignmentType)> {
    stream.clear();

    while let Some((multiplicity, alignment_type)) = alignment.first().copied() {
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            break;
        } else if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
            panic!("Found template switch exit without matching entrance");
        } else {
            stream.push_until_full(multiplicity, alignment_type);
            *alignment = &alignment[1..];

            if stream.is_full() {
                break;
            }
        }
    }

    let downstream = stream.stream_vec();
    stream.pop_to_requested_length();
    downstream
}
