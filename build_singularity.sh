#!/bin/bash
target=${1:-dev}
docker build --target ${target} --tag vast-combined:${target} .
singularity pull --force vast-combined_${target}.sif docker-daemon:vast-combined:${target}
