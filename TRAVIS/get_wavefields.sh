#!/bin/bash

if [ -f "wavefield/bwd_merged/merged_instaseis_db.nc4" ]; then
  echo "Wavefields available from cache"
else
  echo "Wavefields need to be downloaded from Google Drive"

  # 1. Download gdrive client
  wget "https://docs.google.com/uc?id=0B3X9GlR6EmbnQ0FtZmJJUXEyRTA&export=download" -O gdrive
  chmod +x ./gdrive

  # 2. Decrypt gdrive token
  mkdir $HOME/.gdrive
  openssl aes-256-cbc -K $encrypted_50ebf69bd92e_key -iv $encrypted_50ebf69bd92e_iv -in TRAVIS/token_v2.json.enc -out $HOME/.gdrive/token_v2.json -d
  
  # 3. Download wavefield archive with gdrive client
  ./gdrive download 0B2O7p0iBm8A9SWk1WjJZUDR5dVU

  # 4. Untar wavefield files
  tar -xvf kerner_wavefields.tar.bz2

fi
