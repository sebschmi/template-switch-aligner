#!/bin/bash

FILE=$(basename "$1")
DIR=$(dirname "$1")
PREFIX="$DIR/${FILE%.*}"

cargo run --release -- align -p ${PREFIX}.fa -c test_files/config/range -o ${PREFIX}.toml
cargo run --release -- align -p ${PREFIX}.fa -c test_files/config/range -o ${PREFIX}_no_ts.toml --no-ts
cargo run --release -- show --input ${PREFIX}.toml -n ${PREFIX}_no_ts.toml -ps ${PREFIX}.svg