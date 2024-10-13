#!/bin/bash

if ! command -v sage &> /dev/null
then
    exit 77 # no sagemath installation, skip
fi

sage sage-test.sage
