#!/bin/bash
cd ../test-data/TILES/STOKESI_IMAGES
echo Downloading image data from VAST data server.
rsync ubuntu@data.vast-survey.org:/mnt/data/VAST/release/EPOCH38/TILES/STOKESI_IMAGES/image.i.VAST_0334-37.SB50801.cont.taylor.0.restored.conv.fits . -h --progress
rsync ubuntu@data.vast-survey.org:/mnt/data/VAST/release/EPOCH38/TILES/STOKESI_IMAGES/image.i.VAST_0345-31.SB50802.cont.taylor.0.restored.conv.fits . -h --progress
cd ../STOKESI_WEIGHTS
echo Downloading corresponding weight data from VAST data server.
rsync ubuntu@data.vast-survey.org:/mnt/data/VAST/release/EPOCH38/TILES/STOKESI_WEIGHTS/weights.i.VAST_0334-37.SB50801.cont.taylor.0.fits . -h --progress
echo Done. Check your test-data directory. 