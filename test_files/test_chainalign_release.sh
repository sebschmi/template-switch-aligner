#!/bin/bash

cargo run --release -- align --alignment-method a-star-chain-ts -c test_files/config/chainalign --alphabet dna -p "$@"
