#!/bin/bash
if [ -f "wavefield/bwd_merged/merged_instaseis_db.nc4" ]; then
  echo "Wavefields available from cache"
else
  echo "Wavefields need to be downloaded"

  # Check, whether encrypted token for Google Drive account is available
  # Local gdrive installation needs to be activated with it, but since it
  # offers read/write access, it should not be shared with other accounts 
  # than github/seismology.
  # Download from gdrive is much faster (50-100MB/s) than from LMU (<10 MB/s),
  # but pull requests from other github accounts than seismology will not have
  # access to it. 
  if [ $encrypted_50ebf69bd92e_key == "" ]; then

    echo "Downloading from LMU"
    # 1. Download archive from LMU server
    wget https://www.geophysik.uni-muenchen.de/~staehler/kerner_wavefields.tar.bz2

    # 2. Unpack archive
    tar -xvf kerner_wavefields.tar.bz2

  else

    # Download wavefield folder with gdrive client
    ./gdrive download 0BwU9d5SH6pPQNGpkcW9hWHlpYTA --recursive

  fi
fi
