# Check positions of each pol protein
mkdir -p results/12_drm/
python3 sierra-local/scripts/retrieve_hivdb_data.py -size 1000 PR results/12_drm/PR_test.fasta
# PR from 1 to 297
python3 sierra-local/scripts/retrieve_hivdb_data.py -size 1000 RT results/12_drm/RT_test.fasta
# RT goes from 298 to end of aln (698 positions)
# Align seqs to HXB2 for DRM analysis
mafft --6merpair --thread 4 --addfragments data/subtype_seqs_naive/B_UK_final_aln_outliers_removed_final.fasta data/refs/HIV1_REF_2020_pol_DNA_SUBTYPE_B.fasta > results/12_drm/aln_vs_refB_drm.fasta
# Keep 12 positions in the beginning, curate aln, name refB as HXB2 and remove refC
# Break results/12_drm/aln_vs_refB_drm.fasta manually into 2 alns corresponding to PR and RT coordinates
# Replace ? with Ns
awk -F, -v OFS=, '{gsub("?", "N"); print}' results/12_drm/aln_vs_refB_drm_PR.fasta > results/12_drm/aln_vs_refB_drm_PR_adj.fasta
awk -F, -v OFS=, '{gsub("?", "N"); print}' results/12_drm/aln_vs_refB_drm_RT.fasta > results/12_drm/aln_vs_refB_drm_RT_adj.fasta

# Copied env and sierra-local repo from previous run because other ways were not working and usual pip install gets the db v8.8 outdated version
sudo apt install minimap2

# Upgraded to python 3.10
# upgrade pip
curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10

# virtualenv sierra-env 
source sierra-env/bin/activate
# https://github.com/PoonLab/sierra-local/issues/84 (comment some lines)
sudo python3 setup.py install

cd ..
# HIVdb version 9.4 with alignment by post-align
sierralocal results/12_drm/aln_vs_refB_drm_PR_adj.fasta -o results/12_drm/drms_b_PR_sierra.json
sierralocal results/12_drm/aln_vs_refB_drm_RT_adj.fasta -o results/12_drm/drms_b_RT_sierra.json

# Convert json to csv
python3 sierra-local/scripts/json2csv.py results/12_drm/drms_b_PR_sierra.json results/12_drm/drms_b_PR_sierra.csv
python3 sierra-local/scripts/json2csv.py results/12_drm/drms_b_RT_sierra.json results/12_drm/drms_b_RT_sierra.csv

# HTML viewer
# error
# python3 sierra-local/htmlview/HTMLMaker.py results/12_drm/drms_b_PR_sierra.json
# mv results.html results/12_drm/drms_b_PR_sierra.html
# python3 scripts-sierra-local/htmlview/HTMLMaker.py results/12_drm/drms_b_RT_sierra.json 