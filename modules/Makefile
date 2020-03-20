clean:
	rm -rf *.md *.html *.pdf cache/ figure/

%.html: %.Rmd
	./make_slides $(basename $(@))

all: clean 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15



1: 1_nimble_background.html
2: 2_bugs_language.html
3: 3_build_model.html
4: 4_run_mcmc.html
5: 5_customizing_mcmc.html
6: 6_user_dist.html
7: 7_exercises.html
8: 8_meta_analysis.html
9: 9_survival_analysis.html
10: 10_reversible_jump.html
11: 11_operating_model.html
12: 12_model_structure.html
13: 13_nimbleFunctions.html
14: 14_user_sampler.html
15: 15_wrapup.html

