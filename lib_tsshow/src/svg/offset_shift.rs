use lib_tsalign::a_star_aligner::template_switch_distance::TemplateSwitchPrimary;

#[derive(Debug, Default)]
pub struct OffsetShift {
    pub reference: isize,
    pub query: isize,
    /// Translates the index of a source sequence character into the index of a rendered character, accounting for inserted copies before this character.
    reference_source_character_to_rendered_character_index: Vec<usize>,
    /// Translates the index of a source sequence character into the index of a rendered character, accounting for inserted copies before this character.
    query_source_character_to_rendered_character_index: Vec<usize>,
    reference_copies: Vec<TsCopy>,
    query_copies: Vec<TsCopy>,
    current_reference_copy_total_length: usize,
    current_query_copy_total_length: usize,
}

#[derive(Debug)]
pub struct TsCopy {
    remaining_length: usize,
}

pub enum CharacterIsCopy {
    Yes { depth: usize },
    No,
}

impl OffsetShift {
    pub fn push_reference_copy(&mut self, length: usize) {
        assert!(length > 0);
        self.reference_copies.push(TsCopy {
            remaining_length: length,
        });
        self.current_reference_copy_total_length += length;
    }

    pub fn push_query_copy(&mut self, length: usize) {
        assert!(length > 0);
        self.query_copies.push(TsCopy {
            remaining_length: length,
        });
        self.current_query_copy_total_length += length;
    }

    pub fn next_primary(&mut self, primary: TemplateSwitchPrimary) -> CharacterIsCopy {
        let (source_character_to_rendered_character_index, copies, current_copy_total_length) =
            match primary {
                TemplateSwitchPrimary::Reference => (
                    &mut self.reference_source_character_to_rendered_character_index,
                    &mut self.reference_copies,
                    &mut self.current_reference_copy_total_length,
                ),
                TemplateSwitchPrimary::Query => (
                    &mut self.query_source_character_to_rendered_character_index,
                    &mut self.query_copies,
                    &mut self.current_query_copy_total_length,
                ),
            };

        if let Some(TsCopy { remaining_length }) = copies.last_mut() {
            if *remaining_length == 0 {
                copies.pop();
                return self.next_primary(primary);
            }

            *remaining_length -= 1;

            CharacterIsCopy::Yes {
                depth: copies.len().checked_sub(1).unwrap(),
            }
        } else {
            source_character_to_rendered_character_index.push(
                source_character_to_rendered_character_index
                    .last()
                    .copied()
                    .unwrap_or(0)
                    + std::mem::take(current_copy_total_length)
                    + 1,
            );
            CharacterIsCopy::No
        }
    }

    pub fn next_anti_primary(&mut self, primary: TemplateSwitchPrimary) -> CharacterIsCopy {
        self.next_primary(primary.inverted())
    }

    pub fn next_reference(&mut self) -> CharacterIsCopy {
        self.next_primary(TemplateSwitchPrimary::Reference)
    }

    pub fn next_query(&mut self) -> CharacterIsCopy {
        self.next_primary(TemplateSwitchPrimary::Query)
    }

    pub fn primary_source_character_to_rendered_character_index(
        &self,
        primary: TemplateSwitchPrimary,
        source_index: usize,
    ) -> usize {
        match primary {
            TemplateSwitchPrimary::Reference => {
                self.reference_source_character_to_rendered_character_index[source_index]
            }
            TemplateSwitchPrimary::Query => {
                self.query_source_character_to_rendered_character_index[source_index]
            }
        }
    }

    pub fn anti_primary_source_character_to_rendered_character_index(
        &self,
        primary: TemplateSwitchPrimary,
        source_index: usize,
    ) -> usize {
        self.primary_source_character_to_rendered_character_index(primary.inverted(), source_index)
    }
}

impl std::fmt::Display for OffsetShift {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RQ {}/{}", self.reference, self.query)
    }
}
