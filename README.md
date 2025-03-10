# Template Switch Aligner

<table border="0">
 <tr>
    <td>Binary:</td>
    <td><a href="https://crates.io/crates/tsalign"><img src="https://img.shields.io/crates/v/tsalign.svg?style=flat-square" alt="Binary version" /></a></td>
    <td><a href="https://crates.io/crates/tsalign"><img src="https://img.shields.io/crates/d/tsalign.svg?style=flat-square" alt="Binary downloads" /></a></td>
    <td></td>
 </tr>
 <tr>
    <td>Library:</td>
    <td><a href="https://crates.io/crates/lib_tsalign"><img src="https://img.shields.io/crates/v/lib_tsalign.svg?style=flat-square" alt="Library version" /></a></td>
    <td><a href="https://crates.io/crates/lib_tsalign"><img src="https://img.shields.io/crates/d/lib_tsalign.svg?style=flat-square" alt="Library downloads" /></a></td>
    <td><a href="https://docs.rs/lib_tsalign"><img src="https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square" alt="Library docs" /></a></td>
 </tr>
</table>

Align two genomic sequences while allowing for template switches.

## Installation

### Via Cargo (Preferred)

1. Install the rust toolchain by going to [rustup.rs](https://rustup.rs/) and following the instructions.
   Don't worry, on unix-like systems there is just a single command to execute.

2. Run `cargo install tsalign`.

3. You can now run `tsalign` on your command line from anywhere.

If you ever want to update to a new release, simply run `cargo install tsalign` again.

### From Source (For Developers)

1. Install the rust toolchain by going to [rustup.rs](https://rustup.rs/) and following the instructions.
   Don't worry, on unix-like systems there will just be a single command to be executed.

2. Clone this git repository using `git clone <url of this repository>`.

3. From within the root of the git repository, you can run `cargo run --release -- align <arguments>` to run tsalign, where `<arguments>` are the arguments that are passed to tsalign.

Cargo acts as a wrapper here, ensuring that whenever you make changes to the code, it will be recompiled if necessary.
Hence, for updating, it is enough to do a `git pull`.

## Usage

Run the installed tool with `--help` (e.g. `tsalign --help` if installed via cargo) to get an overview of the available options.
