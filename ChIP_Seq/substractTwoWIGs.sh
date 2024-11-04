#!/bin/bash

module add python/3.7.0

python3 /Shared/NEURO/AbelLab/Yann/ChIP-seq/scripts/subtractTwoWig.py /Shared/NEURO/AbelLab/Yann/ChIP-seq/ChIP-seq_PCr_HC_SOR/PCr1.scl.wig /Shared/NEURO/AbelLab/Yann/ChIP-seq/ChIP-seq_PCr_HC_SOR/Input1.scl.wig /Shared/NEURO/AbelLab/Yann/ChIP-seq/ChIP-seq_PCr_HC_SOR/PCr1.sclCtl.wig
