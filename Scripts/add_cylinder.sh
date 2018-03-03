#!/bin/bash
imp=$HOME/imp_git/fast/setup_environment.sh
imppy="$imp /usr/bin/python"
if [ -e $2 ]; then
    echo "$2 already exists - as a safety measure, not overwriting. Delete explicitly if needed, after making sure you didn't mix input $1 with output $2"
    exit -1
fi
$imppy ~/npctransport/utility/add_cylinder.py --input_rmf $1 --output_rmf $2 --ref_output $3 --smooth_n_frames=100 --radius=2 --recolor_fgs --skip_n_frames=50
