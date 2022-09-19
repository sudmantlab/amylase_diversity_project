
make 3 files
--ancientSamps [time series of genotype likelihoods]
    - age (in generations), likelihood
--coal inferred european population size
--popFreq single epoch spanning 15000 years

~/Documents/science/programs/clues/inference.py


coal/coalescence_rates are from relate paper


python inference.py --times example/example --ancientSamps example/exampleAncientSamps.txt

exampleAncientSamps.txt is a file listed ancient genotype likelihoods (set likelihoods to 0/-inf if hard-called)


/Users/petersudmant/Documents/science/sudmantlab/projects/ancient_DNA/amylase/amylase_diversity_project/haplotype_analysis/clues/output/ancient_timecourse/AMY/data/104097300.AncientSamps.txt



ADD IN MODERN FREQ - SEE WHAT HAPPENS


