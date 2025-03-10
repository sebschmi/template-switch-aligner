#![deny(clippy::mod_module_files)]

use clap::Parser;

mod align;

#[derive(Parser)]

struct Cli {
    #[command(subcommand)]
    subcommand: Subcommand,
}

#[derive(clap::Subcommand)]
enum Subcommand {
    Align(align::Cli),
}

fn main() {
    let cli = Cli::parse();

    match cli.subcommand {
        Subcommand::Align(cli) => align::align_cli(cli),
    }
}
