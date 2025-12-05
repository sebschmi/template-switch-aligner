#!/bin/bash

cargo run -- align --alignment-method a-star-chain-ts -c test_files/config/chainalign --alphabet dna -p "$@"