#!/bin/sh

# Script to update local versions of the metadata

echo "\nMETADATA"
gdown 'https://drive.google.com/uc?export=download&id=1qOzys6wmWXGDDwgGs46XF4GGNarq3_JG'
mv _metadata_saccharomyces_cerevisiae.xlsx ~/winhome/Desktop/Projets/sourdough_evolution/data/_metadata_saccharomyces_cerevisiae.xlsx

echo "\nSTRAINS"
gdown 'https://drive.google.com/uc?export=download&id=1KCI63qEKIZtpADd9bRkHEi_vC1766TW4'
mv strains.xlsx ~/winhome/Desktop/Projets/sourdough_evolution/data/strains.xlsx

echo "\nBAKERS"
gdown 'https://drive.google.com/uc?export=download&id=1Uq5o-wpCG_IcTDcxoE6sblQJCEbLKahT'
mv bakers.xlsx ~/winhome/Desktop/Projets/sourdough_evolution/data/bakers.xlsx

echo "\nFLOURS"
gdown 'https://drive.google.com/uc?export=download&id=1EfZ-MjVZnXt2hdmegpVIOdhY5Ywz_VGI'
mv flours.xlsx ~/winhome/Desktop/Projets/sourdough_evolution/data/flours.xlsx

echo "\nDATA CYTO"
gdown 'https://drive.google.com/uc?export=download&id=14SRs-YB4zQPcSDkLS4Gf1Z7_enXEGhVi'
mv data_cyto.xlsx ~/winhome/Desktop/Projets/sourdough_evolution/data/data_robot/data_cyto.xlsx

echo "\nDATA OTHER_STRAINS"
gdown 'https://drive.google.com/uc?export=download&id=1UsCiDjrLvbVSmW3xb8LHcrrjqpPv2twI'
mv other_strains_metadata.xlsx ~/winhome/Desktop/Projets/sourdough_evolution/data/other_strains/other_strains_metadata.xlsx
