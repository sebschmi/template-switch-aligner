use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::GenomeSequence,
};
use lib_ts_chainalign::{
    chain_align::AlignmentParameters, chaining_lower_bounds::ChainingLowerBounds,
    costs::AlignmentCosts,
};
use lib_tsalign::{
    a_star_aligner::alignment_geometry::AlignmentRange, config::TemplateSwitchConfig,
};
use log::{debug, info, warn};
use sha::{
    sha1::Sha1,
    utils::{Digest, DigestExt},
};
use std::{fmt::Debug, fs::File, path::PathBuf};

use crate::align::Cli;

pub fn align_a_star_chain_ts<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
) {
    let mut config_path = cli.configuration_directory.clone();
    info!("Loading alignment config directory {config_path:?}");

    config_path.push("config.tsa");
    let config_file = std::io::BufReader::new(
        std::fs::File::open(&config_path)
            .unwrap_or_else(|error| panic!("Error opening config file {config_path:?}: {error}")),
    );
    let alignment_costs = TemplateSwitchConfig::<AlphabetType, _>::read_plain(config_file)
        .unwrap_or_else(|error| panic!("Error parsing template switch config:\n{error}"));
    let alignment_costs: AlignmentCosts<_> = alignment_costs.into();

    let cache_directory = cli.cache_directory.clone().unwrap_or_else(|| {
        warn!("No cache directory specified, dropping files into current working directory.");
        PathBuf::new()
    });
    let max_n = 1 << (usize::BITS - (reference.len().max(query.len()) - 1).leading_zeros());
    let k = cli.k.unwrap_or_else(|| {
        // This evaluates to ceil(log_2(length_sum)).
        // The motivation is that there are length_sum k-mers,
        // so for each to be different, k needs to be at least ceil(log_4(length_sum)).
        // However, the birthday paradoxon states that for avoiding collisions,
        // the amount of possible k-mers needs to grow in the square of the amount available k-mers,
        // so we square that and arrive at ceil(log_2(length_sum)).
        let k = usize::BITS - ((reference.len() + query.len()) - 1).leading_zeros();
        // Decrease k a little, because we can hopefully afford a few more anchors.
        k.saturating_sub(3).max(2)
    });
    debug!("Using max_n = {max_n}");
    info!("Using k = {k}");
    let max_match_run = k - 1;
    let cost_hash = Sha1::default()
        .digest(
            &bincode::serde::encode_to_vec(&alignment_costs, bincode::config::standard()).unwrap(),
        )
        .to_hex();
    debug!("Using cost_hash = {cost_hash}");

    let cache_file: PathBuf = [
        cache_directory,
        PathBuf::from(format!("{cost_hash}-{k}-{max_n}.tsc")),
    ]
    .iter()
    .collect();

    let chaining_lower_bounds = if let Ok(mut file) = File::open(&cache_file) {
        info!("Loading preprocessed data from cache at {cache_file:?}");
        let chaining_lower_bounds = ChainingLowerBounds::read(&mut file).unwrap();
        assert_eq!(chaining_lower_bounds.alignment_costs(), &alignment_costs);
        assert_eq!(chaining_lower_bounds.max_match_run(), max_match_run);
        chaining_lower_bounds
    } else {
        info!("Preprocessing...");
        let chaining_lower_bounds =
            lib_ts_chainalign::preprocess(max_n, max_match_run, alignment_costs);

        info!("Storing preprocessed data into cache at {cache_file:?}");
        let mut file = File::create(&cache_file).unwrap();
        chaining_lower_bounds.write(&mut file).unwrap();
        chaining_lower_bounds
    };

    let reference = reference.clone_as_vec();
    let query = query.clone_as_vec();
    let parameters = AlignmentParameters {
        max_successors: cli.max_chaining_successors,
        max_exact_cost_function_cost: cli.max_exact_cost_function_cost.into(),
        closed_list: cli.chaining_closed_list,
        open_list: cli.chaining_open_list,
    };

    let alignment = lib_ts_chainalign::align::<AlphabetType>(
        reference,
        query,
        range,
        &parameters,
        &|c| {
            AlphabetType::character_to_ascii(
                AlphabetType::ascii_to_character(c).unwrap().complement(),
            )
        },
        reference_name,
        query_name,
        &chaining_lower_bounds,
    );
    info!("Finished aligning");

    if let Some(output) = cli.output {
        info!("Outputting alignment statistics to {output:?}");
        use std::io::Write;
        let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        write!(output, "{}", toml::to_string(&alignment).unwrap()).unwrap();
    }

    println!("{alignment}");
}
