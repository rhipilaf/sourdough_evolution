#!/bin/bash

touch check

which truc &> check
OK=$?

rm check

if [ $OK -eq 0 ]; then
  echo "\ngdown needs to be installed to download the file with this command : \n   sudo pip install gdown\n\nThen re-run this script. It should work."
  exit 1
fi

#The following downloads the up-to-date xlsx file that contains metadata from Google drive and puts it in the data/ folder

gdown 'https://drive.google.com/uc?export=download&id=1qOzys6wmWXGDDwgGs46XF4GGNarq3_JG'
mv metadata_saccharomyces_cerevisiae.xlsx ~/winhome/Desktop/Projets/sourdough_evolution/data/metadata_saccharomyces_cerevisiae.xlsx