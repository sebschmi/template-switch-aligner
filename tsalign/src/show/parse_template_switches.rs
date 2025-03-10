use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;

use super::alignment_stream::{AlignmentCoordinates, AlignmentStream};

pub struct TSShow<AlignmentType> {
    pub offset: AlignmentCoordinates,
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
    let offset = stream.tail_coordinates();
    let upstream = stream.stream_vec();
    let mut template_switch = Vec::new();

    while let Some((multiplicity, alignment_type)) = alignment.first().copied() {
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            panic!("Found template switch entrance within template switch");
        } else if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
            template_switch.push((multiplicity, alignment_type));
            stream.push(multiplicity, alignment_type);
            *alignment = &alignment[1..];

            let template_switch = template_switch;
            let downstream = parse_downstream(alignment, stream);

            return TSShow {
                offset,
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
