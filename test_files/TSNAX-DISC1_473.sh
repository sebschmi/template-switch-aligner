#!/usr/bin/env bash

cargo build --release

target/release/tsalign align -p test_files/TSNAX-DISC1_473.fa --alignment-method a-star-template-switch --ts-total-length-strategy none --skip-characters '-' --configuration-directory test_files/config/bench/ --alphabet dna --rq-ranges R196..219Q196..212 --ts-min-length-strategy none
target/release/tsalign align -p test_files/TSNAX-DISC1_473.fa --alignment-method a-star-template-switch --ts-total-length-strategy none --skip-characters '-' --configuration-directory test_files/config/bench/ --alphabet dna --rq-ranges R196..219Q196..212 --ts-min-length-strategy lookahead
target/release/tsalign align -p test_files/TSNAX-DISC1_473.fa --alignment-method a-star-template-switch --ts-total-length-strategy none --skip-characters '-' --configuration-directory test_files/config/bench/ --alphabet dna --rq-ranges R196..219Q196..212 --ts-min-length-strategy preprocess-price
target/release/tsalign align -p test_files/TSNAX-DISC1_473.fa --alignment-method a-star-template-switch --ts-total-length-strategy none --skip-characters '-' --configuration-directory test_files/config/bench/ --alphabet dna --rq-ranges R196..219Q196..212 --ts-min-length-strategy preprocess-lookahead