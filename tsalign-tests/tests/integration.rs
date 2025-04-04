use anyhow::Result;
use util::run_in_repo_root;

mod util;

#[test]
fn test_align_default_cfg_twin() -> Result<()> {
    run_in_repo_root("align -p test_files/twin_a.fa")
}

#[ignore = "Currently doesn't work 🫣"]
#[test]
fn test_align_default_cfg_qr() -> Result<()> {
    run_in_repo_root("align -q test_files/query_a.fa -r test_files/reference_a.fa")
}
