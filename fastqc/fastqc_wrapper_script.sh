#!/bin/bash
i=$1
fastqc ${i}/${i}_R1_trimmed.fq
fastqc ${i}/${i}_R2_trimmed.fq
