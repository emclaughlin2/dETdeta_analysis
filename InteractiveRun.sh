#!/bin/bash

# $1 = runnumber, $2 = generator, $3 = minus_z, $4 = plus_z, $5 = dataormc, $6 = reweighting, $7 = central
root.exe -q -b "dETdeta_vertex_reweighting.C(\"$1\", \"$2\")"
root.exe -q -b "dETdeta_hot_dead_map.C(\"$1\", \"$3\", \"$4\")"
root.exe -q -b "dETdeta_analysis.C(\"$1\", \"$2\", \"$3\", \"$4\", \"$5\", \"$6\", \"$7\")"
echo "JOB COMPLETE!"
