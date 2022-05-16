#!/bin/bash
i=$1
fastqc ${i}/${i}_R1.fq.gz
fastqc ${i}/${i}_R2.fq.gz
