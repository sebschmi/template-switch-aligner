#![deny(clippy::mod_module_files)]

use clap::Parser;

mod align;
mod show;

#[derive(Parser)]

struct Cli {
    #[command(subcommand)]
    subcommand: Subcommand,
}

#[derive(clap::Subcommand)]
enum Subcommand {
    Align(align::Cli),
    Show(show::Cli),
}

fn main() {
    let cli = Cli::parse();

    match cli.subcommand {
        Subcommand::Align(cli) => align::cli(cli),
        Subcommand::Show(cli) => show::cli(cli),
    }
}
