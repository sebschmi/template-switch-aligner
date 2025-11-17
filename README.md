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

## Features

* Align two genomic sequences while allowing for template switches.
* Visualise the alignment with template switch jumps.


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


### Aligning

Tsalign takes as input either a single fasta file with two fasta records, or two separate fasta files.
The output is a toml file which contains alignment statistics and can be used to visualise the alignment.

The alignment metric can be configured via a config file named `config.tsa`.
The path to the file can be passed via the `-c` parameter.
Only pass the directory that contains the file, and not the file itself.

Note that tsalign is very resource-intensive, so it often makes sense to use the `--memory-limit` parameter, which takes a memory limit in bytes.
This limit is applied approximately, with an accuracy of about a factor of two.

If you have two sequences in a file `pair.fa`, computing an alignment could look as follows.

```bash
tsalign align -p pair.fa -o alignment.toml -c sample_tsa_config --memory-limit 1000000000
```

#### Choosing a different alphabet

If you want to align strings with an alphabet other than the default DNA alphabet, consider using the parameter `-a`.
Remember to also update the `config.tsa` file with a metric that corresponds to the characters available in the alphabet.
Available alphabets are:

* **dna**: ACGT
* **dna-n**: ACGTN
* **rna**: ACGU
* **rna-n**: ACGUN
* **dna-iupac**: ABCDGHKMNRSTVWY
* **rna-iupac**: ABCDGHKMNRSUVWY

### Visualisation

### Setting the alignment range

When searching for template switches, researchers are often interested in specific mutation hotspots that are only tens of characters long, while the 2-3-alignment of the template switch may also align outside of the mutation hotspot.
This requires to pass longer sequences to tsalign, to make the context around the mutation hotspot available for alignment.
However, since tsalign's performance is very sensitive to the lengths of the input sequences, it is important to specify the range of the mutation hotspot, such that tsalign computes the alignment only of this shorter region, while only allowing 2-3-alignments of template switches to align outside of the region.

For an example, consider the following input sequences:

Original `input.fa`.
```fasta
>reference
ACACACCCAACGCGGG
>query
ACAAACGTGTCGCGCG
```

We want to check if the substitution of `CCAA` to `GTGT` can be better explained with a template switch.
Therefore, we want tsalign to focus on this region.
We can insert `|` characters to mark the focus region, including an extra character at the beginning and the end for good measure.

Modified `input.delimited.fa`.
```fasta
>reference
ACACA|CCCAAC|GCGGG
>query
ACAAA|CGTGTC|GCGCG
```

Now we can run `tsalign align -p input.delimited.fa --use-embedded-rq-ranges`.
Now, tsalign will use much less resources, as it can ignore the non-matches before and after the focus region.
