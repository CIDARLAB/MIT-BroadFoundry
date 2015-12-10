
BIN_PATH=/Users/Tom/Dropbox/Research/projects/mit_voigt_lab/02_circuit_analysis/_support/goatools/scripts

python make_studies.py

python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/flask_vs_tube.down.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/flask_vs_tube.down.go.txt
python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/flask_vs_tube.up.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/flask_vs_tube.up.go.txt

python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/broken_flask_vs_tube.down.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/broken_flask_vs_tube.down.go.txt
python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/broken_flask_vs_tube.up.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/broken_flask_vs_tube.up.go.txt

python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/ara_comp_tube.down.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/ara_comp_tube.down.go.txt
python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/ara_comp_tube.up.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/ara_comp_tube.up.go.txt

python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/iptg_comp_tube.down.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/iptg_comp_tube.down.go.txt
python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/iptg_comp_tube.up.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/iptg_comp_tube.up.go.txt

python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/atc_comp_tube.down.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/atc_comp_tube.down.go.txt
python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/atc_comp_tube.up.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/atc_comp_tube.up.go.txt

python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/rep_3_yfp_tube.down.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/rep_3_yfp_tube.down.go.txt
python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/rep_3_yfp_tube.up.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/rep_3_yfp_tube.up.go.txt

python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/rep_4_tube.down.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/rep_4_tube.down.go.txt
python $BIN_PATH/find_enrichment.py --pval=0.01 --indent ./studies/rep_4_tube.up.txt ./data/ecoli_go_population.txt ./data/ecoli_go_association.txt > ./results/rep_4_tube.up.go.txt
