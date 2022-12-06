#!/bin/bash -l

release_epoch=${1:-EPOCH00}
data_dir=${2:-neighbours}

rsync -vhsrLt -e ssh --progress \
    --include "${release_epoch}/" \
    --include "${release_epoch}/RACS_*/" \
    --include "${release_epoch}/RACS_*/RACS_*.EPOCH*.[IV].conv.fits" \
    --include "${release_epoch}/RACS_*/RACS_*.EPOCH*.[IV].conv.weight.fits" \
    --include "${release_epoch}/RACS_*/noiseMap.RACS_*.EPOCH*.[IV].conv.fits" \
    --include "${release_epoch}/RACS_*/meanMap.RACS_*.EPOCH*.[IV].conv.fits" \
    --include "${release_epoch}/RACS_*/*.xml" \
    --exclude "*" \
    ${data_dir}/${release_epoch} ubuntu@data.vast-survey.org:/data/.staging/combined/
