#!/bin/bash

DIR=/sf/bd/applications/photon_spectra

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

export PYTHONPATH=$DIR/pythonpath/uncertainties/:$DIR/pythonpath/:/sf/bd/members/Philipp/pythonpath/slic/

cd $DIR/gui
python main.py --facility SwissFEL

