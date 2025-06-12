#!/usr/bin/env bash

cargo run -- align -p test_files/twin_min_length_test.fa --rq-ranges R40..75Q41..75
cargo run -- align -p test_files/twin_min_length_test.fa --rq-ranges R40..75Q41..75 --ts-min-length-strategy preprocess-price
cargo run -- align -p test_files/twin_min_length_test.fa --rq-ranges R40..75Q41..75 --ts-min-length-strategy preprocess-filter