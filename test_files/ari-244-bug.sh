#!/usr/bin/env bash

nohup cargo run -- -p test_files/twin_ari_email_244.fa --skip-characters - --alphabet dna-n -c test_files/config/experiments --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy lookahead >ari-244-bug.log 2>&1 &