use std::{
    fs::File,
    io::{stdout, Read, Write},
    path::PathBuf,
};

use clap::Parser;

use lib_tsshow::{plain_text::show_template_switches, svg::create_ts_svg, svg_to_png};
use log::{info, LevelFilter};
use simplelog::{ColorChoice, TermLogger, TerminalMode};

#[derive(Parser)]
pub struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    /// Path to a toml output file of `tsalign align`.
    ///
    /// This is the input from which template switches are visualised.
    #[clap(long, short = 'i')]
    input: PathBuf,

    /// Path to a toml output file of `tsalign align` created with the `--no-ts` flag.
    ///
    /// This input can be used to provide an alternative alignment without template switches.
    /// It is expected to contain no template switches, otherwise it will be rejected.
    #[clap(long, short = 'n')]
    no_ts_input: Option<PathBuf>,

    /// Create an SVG image file showing the complete alignment with all template switches.
    #[clap(long, short = 's')]
    svg: Option<PathBuf>,

    /// If set together with --svg, then the SVG image is rendered to a PNG image as well.
    #[clap(long, short = 'p')]
    png: bool,
}

pub fn cli(cli: Cli) {
    TermLogger::init(
        cli.log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )
    .unwrap();

    info!("Reading tsalign output toml file {:?}", cli.input);
    let mut buffer = String::new();
    File::open(cli.input)
        .unwrap_or_else(|error| panic!("Error opening input file: {error}"))
        .read_to_string(&mut buffer)
        .unwrap_or_else(|error| panic!("Error reading input file: {error}"));
    let result =
        toml::from_str(&buffer).unwrap_or_else(|error| panic!("Error parsing input file: {error}"));

    let no_ts_result = cli.no_ts_input.as_ref().map(|no_ts_input| {
        info!("Reading tsalign no-ts output toml file {no_ts_input:?}");
        buffer.clear();
        File::open(no_ts_input)
            .unwrap_or_else(|error| panic!("Error opening input file: {error}"))
            .read_to_string(&mut buffer)
            .unwrap_or_else(|error| panic!("Error reading input file: {error}"));
        toml::from_str(&buffer).unwrap_or_else(|error| panic!("Error parsing input file: {error}"))
    });

    show_template_switches(stdout(), &result, &no_ts_result);

    if let Some(svg_out_path) = cli.svg.as_ref() {
        let mut svg = Vec::new();
        create_ts_svg(&mut svg, &result, &no_ts_result);
        let svg = svg;

        info!("Writing svg to {svg_out_path:?}");
        File::create(svg_out_path).unwrap().write_all(&svg).unwrap();

        if cli.png {
            let png = svg_to_png(&svg, 20.0);

            let png_out_path = svg_out_path.with_extension("png");
            info!("Writing png to {png_out_path:?}");
            File::create(png_out_path).unwrap().write_all(&png).unwrap();
        }
    }
}
