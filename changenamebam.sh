#!/bin/bash

samtools view -H $1 | sed  -e 's/SN:\([0-9XY]*\)/SN:judiseq_\1/' -e 's/SN:MT/SN:judiseq_M/' | samtools reheader - $1 > $2
