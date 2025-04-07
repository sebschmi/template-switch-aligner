use std::{
    env,
    process::{Command, Stdio},
};

use anyhow::{Result, anyhow};

pub fn run_in_repo_root(args: &str) -> Result<()> {
    if env::current_dir()?.ends_with("tsalign-tests") {
        // working directory is this crate, a.k.a. "[...]/template-switch-aligner/tsalign-tests"
        // simulate a call from the repo root by traversing to "../"
        env::set_current_dir(
            env::current_dir()?
                .parent()
                .ok_or(anyhow!("No parent directory"))?,
        )?;
    }

    let process = Command::new("cargo")
        .args(["run", "--"].into_iter().chain(args.split_whitespace()))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    let output = process.wait_with_output()?;

    println!("{}", output.status);
    println!("stdout:\n{}", String::from_utf8_lossy(&output.stdout));
    println!("stderr:\n{}", String::from_utf8_lossy(&output.stderr));

    assert!(output.status.success());

    Ok(())
}
