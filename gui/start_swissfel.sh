#!/bin/bash

source /opt/gfa/python 3.7


__conda_setup="$('/sf/bd/packages/conda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/sf/bd/packages/conda/etc/profile.d/conda.sh" ]; then
        . "/sf/bd/packages/conda/etc/profile.d/conda.sh"
    else
        export PATH="/sf/bd/packages/conda/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate slic
export PYTHONPATH=../pythonpath/uncertainties/:../pythonpath/:/sf/bd/members/Philipp/pythonpath/slic/

python main.py --facility SwissFEL

