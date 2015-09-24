
python ../../_support/goatools/scripts/find_enrichment.py --pval=0.01 --indent go_results/input_678_down.txt ecoli_go_population.txt ecoli_go_association.txt > go_results/678_down.results.txt

python ../../_support/goatools/scripts/find_enrichment.py --pval=0.01 --indent go_results/input_678_up.txt ecoli_go_population.txt ecoli_go_association.txt > go_results/678_up.results.txt
