#!/usr/bin/env bash

cargo build

if [ ! -f test_files/STXBP5-AS1_99.toml ]; then
    cargo build --release
    target/release/tsalign align -p test_files/STXBP5-AS1_99.fa -o test_files/STXBP5-AS1_99.toml --alignment-method a-star-template-switch --skip-characters 'N-' --alphabet dna -c test_files/config/bench --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy preprocess-price -k 0 --max-chaining-successors 0 --max-exact-cost-function-cost 0 --chaining-closed-list special --chaining-open-list linear-heap --rq-ranges R196..227Q196..202
    target/release/tsalign align -p test_files/STXBP5-AS1_99.fa -o test_files/STXBP5-AS1_99.no_ts.toml --alignment-method a-star-template-switch --skip-characters 'N-' --alphabet dna -c test_files/config/bench --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy preprocess-price -k 0 --max-chaining-successors 0 --max-exact-cost-function-cost 0 --chaining-closed-list special --chaining-open-list linear-heap --rq-ranges R196..227Q196..202 --no-ts
fi

target/debug/tsalign show -i test_files/STXBP5-AS1_99.toml -n test_files/STXBP5-AS1_99.no_ts.toml -ps test_files/STXBP5-AS1_99.svg --log-level trace
