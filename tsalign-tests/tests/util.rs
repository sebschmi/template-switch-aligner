use std::env;

use anyhow::{Result, anyhow};
use clap::Parser;
use tsalign::{align, show};

pub fn run_in_repo_root(args: &str) -> Result<()> {
    // working directory is this crate, a.k.a. "[...]/template-switch-aligner/tsalign-tests"
    // simulate a call from the repo root by traversing to "../"
    env::set_current_dir(
        env::current_dir()?
            .parent()
            .ok_or(anyhow!("No parent directory"))?,
    )?;

    if args.starts_with("align ") {
        let args = align::Cli::parse_from(args.split_whitespace());
        align::cli(args)?;
    } else if args.starts_with("show ") {
        let args = show::Cli::parse_from(args.split_whitespace());
        show::cli(args)?;
    }

    Ok(())
}
