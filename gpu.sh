#!/bin/bash

if command -v nvidia-smi &> /dev/null
then
    echo "NVIDIA GPU detected. Current status:"
    nvidia-smi
else
    echo "NVIDIA GPU or nvidia-smi not found."