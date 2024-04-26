#!/bin/bash

for i in {0..123}
do
    ./analysis ${i}
    echo "Job ${i} done!"
done