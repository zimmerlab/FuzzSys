PYTHON = python3
FUZZSYS = PATH_TO_FUZZSYS_DIRECTORY


test1_dataset: test1_dataset_prepare test1_dataset_fuzzy test1_dataset_report
test1_normalweight: test1_normalweight_prepare test1_normalweight_fuzzy test1_normalweight_report
test1_overweight: test1_overweight_prepare test1_overweight_fuzzy test1_overweight_report
test1: test1_dataset test1_normalweight test1_overweight
test2_dataset: test2_dataset_prepare test2_dataset_fuzzy test2_dataset_report
test2_adult: test2_adult_prepare test2_adult_fuzzy test2_adult_report
test2_kid: test2_kid_prepare test2_kid_fuzzy test2_kid_report
test2: test2_dataset test2_adult test2_kid test2_joint_report
test3: test3_prepare test3_fuzzy test3_report
test4: test4_prepare test4_fuzzy test4_report
test5: test5_prepare test5_fuzzy test5_report
test6: test6_prepare test6_fuzzy test6_report
test7_dataset: test7_dataset_prepare test7_dataset_fuzzy test7_dataset_report
test7_normalweight: test7_normalweight_prepare test7_normalweight_fuzzy test7_normalweight_report
test7_overweight: test7_overweight_prepare test7_overweight_fuzzy test7_overweight_report
test7: test7_dataset test7_normalweight test7_overweight test7_joint_report
run_test: test1 test2 test3 test4 test5 test6 test7
.PHONY: prepare test1 test2 test3 test4 test5 test6 test7


prepare:	prepare_input.py \
			$(wildcard ./excel_tables/*.xlsx)
	$(PYTHON) prepare_input.py --excel ./excel_tables/ \
		--base 2 \
		--pseudo 0 \
		--output ./data/


test1_dataset_prepare:	prepare_submatrix.py \
						./data/raw_expression.tsv \
						./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "underweight,normalweight,overweight,class I obesity,class II obesity" \
		--output ./test1/data_set/adult.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-N,kUmNA-N,kNmA-N,kNmNA-N" \
		--output ./test1/data_set/kid.tsv


test1_dataset_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						./test1/data_set/adult.tsv \
						./test1/data_set/kid.tsv \
						config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test1/data_set/adult.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test1/data_set/fuzzy_adult/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test1/data_set/kid.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test1/data_set/fuzzy_kid/


test1_dataset_report:	generate_report.py \
						./test1/data_set/adult.tsv \
						./test1/data_set/kid.tsv \
						./test1/data_set/fuzzy_adult/concepts_detailed.json \
						./test1/data_set/fuzzy_kid/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test1/data_set/ \
		--label1 adult \
		--label2 kid \
		--output ./test1/data_set/comparison_reports/


test1_normalweight_prepare:	prepare_submatrix.py \
							./data/raw_expression.tsv \
							./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "underweight,normalweight" \
		--output ./test1/normalweight/adult.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kNmA-N,kNmNA-N" \
		--output ./test1/normalweight/kid.tsv


test1_normalweight_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
							./test1/normalweight/adult.tsv \
							./test1/normalweight/kid.tsv \
							config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test1/normalweight/adult.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test1/normalweight/fuzzy_adult/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test1/normalweight/kid.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test1/normalweight/fuzzy_kid/


test1_normalweight_report:	generate_report.py \
							./test1/normalweight/adult.tsv \
							./test1/normalweight/kid.tsv \
							./test1/normalweight/fuzzy_adult/concepts_detailed.json \
							./test1/normalweight/fuzzy_kid/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test1/normalweight/ \
		--label1 adult \
		--label2 kid \
		--output ./test1/normalweight/comparison_reports/


test1_overweight_prepare:	prepare_submatrix.py \
							./data/raw_expression.tsv \
							./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "overweight,class I obesity,class II obesity" \
		--output ./test1/overweight/adult.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-N,kUmNA-N" \
		--output ./test1/overweight/kid.tsv


test1_overweight_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						./test1/overweight/adult.tsv \
						./test1/overweight/kid.tsv \
						config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test1/overweight/adult.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test1/overweight/fuzzy_adult/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test1/overweight/kid.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test1/overweight/fuzzy_kid/


test1_overweight_report:	generate_report.py \
							./test1/overweight/adult.tsv \
							./test1/overweight/kid.tsv \
							./test1/overweight/fuzzy_adult/concepts_detailed.json \
							./test1/overweight/fuzzy_kid/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test1/overweight/ \
		--label1 adult \
		--label2 kid \
		--output ./test1/overweight/comparison_reports/


test2_dataset_prepare:	prepare_submatrix.py \
						./data/raw_expression.tsv \
						./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "overweight,class I obesity,class II obesity,kUmA-N,kUmNA-N" \
		--output ./test2/data_set/overweight.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "underweight,normalweight,kNmA-N,kNmNA-N" \
		--output ./test2/data_set/normalweight.tsv


test2_dataset_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						./test2/data_set/overweight.tsv \
						./test2/data_set/normalweight.tsv \
						config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test2/data_set/overweight.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test2/data_set/fuzzy_overweight/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test2/data_set/normalweight.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test2/data_set/fuzzy_normalweight/


test2_dataset_report:	generate_report.py \
						./test2/data_set/overweight.tsv \
						./test2/data_set/normalweight.tsv \
						./test2/data_set/fuzzy_overweight/concepts_detailed.json \
						./test2/data_set/fuzzy_normalweight/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test2/data_set/ \
		--label1 overweight \
		--label2 normalweight \
		--output ./test2/data_set/comparison_reports/


test2_adult_prepare:	prepare_submatrix.py \
						./data/raw_expression.tsv \
						./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "overweight,class I obesity,class II obesity" \
		--output ./test2/adult/overweight.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "underweight,normalweight" \
		--output ./test2/adult/normalweight.tsv


test2_adult_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
					./test2/adult/overweight.tsv \
					./test2/adult/normalweight.tsv \
					config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test2/adult/overweight.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test2/adult/fuzzy_overweight/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test2/adult/normalweight.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test2/adult/fuzzy_normalweight/


test2_adult_report:	generate_report.py \
					./test2/adult/overweight.tsv \
					./test2/adult/normalweight.tsv \
					./test2/adult/fuzzy_overweight/concepts_detailed.json \
					./test2/adult/fuzzy_normalweight/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test2/adult/ \
		--label1 overweight \
		--label2 normalweight \
		--output ./test2/adult/comparison_reports/


test2_kid_prepare:	prepare_submatrix.py \
					./data/raw_expression.tsv \
					./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-N,kUmNA-N" \
		--output ./test2/kid/overweight.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kNmA-N,kNmNA-N" \
		--output ./test2/kid/normalweight.tsv


test2_kid_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
					./test2/kid/overweight.tsv \
					./test2/kid/normalweight.tsv \
					config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test2/kid/overweight.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test2/kid/fuzzy_overweight/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test2/kid/normalweight.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test2/kid/fuzzy_normalweight/


test2_kid_report:	generate_report.py \
					./test2/kid/overweight.tsv \
					./test2/kid/normalweight.tsv \
					./test2/kid/fuzzy_overweight/concepts_detailed.json \
					./test2/kid/fuzzy_normalweight/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test2/kid/ \
		--label1 overweight \
		--label2 normalweight \
		--output ./test2/kid/comparison_reports/


test2_joint_report:	generate_report.py \
					./test2/adult/overweight.tsv \
					./test2/adult/normalweight.tsv \
					./test2/kid/overweight.tsv \
					./test2/kid/normalweight.tsv \
					./test2/adult/fuzzy_overweight/concepts_detailed.json \
					./test2/adult/fuzzy_normalweight/concepts_detailed.json \
					./test2/kid/fuzzy_overweight/concepts_detailed.json \
					./test2/kid/fuzzy_normalweight/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test2/ \
		--label1 overweight \
		--label2 normalweight \
		--group1 adult \
		--group2 kid \
		--output ./test2/joint_reports/


test3_prepare:	prepare_submatrix.py \
				./data/raw_expression.tsv \
				./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-N" \
		--output ./test3/kUmA-N.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmNA-N,kNmA-N,kNmNA-N" \
		--output ./test3/kXmX-N.tsv


test3_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
				./test3/kUmA-N.tsv \
				./test3/kXmX-N.tsv \
				config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test3/kUmA-N.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test3/fuzzy_kUmA-N/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test3/kXmX-N.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test3/fuzzy_kXmX-N/


test3_report:	generate_report.py \
				./test3/kUmA-N.tsv \
				./test3/kXmX-N.tsv \
				./test3/fuzzy_kUmA-N/concepts_detailed.json \
				./test3/fuzzy_kXmX-N/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test3/ \
		--label1 kUmA-N \
		--label2 kXmX-N \
		--output ./test3/comparison_reports/


test4_prepare:	prepare_submatrix.py \
				./data/raw_expression.tsv \
				./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-N,kUmNA-N" \
		--output ./test4/kUmX-N.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kNmA-N,kNmNA-N" \
		--output ./test4/kNmX-N.tsv


test4_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
				./test4/kUmX-N.tsv \
				./test4/kNmX-N.tsv \
				config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test4/kUmX-N.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test4/fuzzy_kUmX-N/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test4/kNmX-N.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test4/fuzzy_kNmX-N/


test4_report:	generate_report.py \
				./test4/kUmX-N.tsv \
				./test4/kNmX-N.tsv \
				./test4/fuzzy_kUmX-N/concepts_detailed.json \
				./test4/fuzzy_kNmX-N/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test4/ \
		--label1 kUmX-N \
		--label2 kNmX-N \
		--output ./test4/comparison_reports/


test5_prepare:	prepare_submatrix.py \
				./data/raw_expression.tsv \
				./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-T" \
		--output ./test5/kUmA-T.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmNA-T,kNmA-T,kNmNA-T" \
		--output ./test5/kXmX-T.tsv


test5_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
				./test5/kUmA-T.tsv \
				./test5/kXmX-T.tsv \
				config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test5/kUmA-T.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test5/fuzzy_kUmA-T/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test5/kXmX-T.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test5/fuzzy_kXmX-T/


test5_report:	generate_report.py \
				./test5/kUmA-T.tsv \
				./test5/kXmX-T.tsv \
				./test5/fuzzy_kUmA-T/concepts_detailed.json \
				./test5/fuzzy_kXmX-T/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test5/ \
		--label1 kUmA-T \
		--label2 kXmX-T \
		--output ./test5/comparison_reports/


test6_prepare:	prepare_submatrix.py \
				./data/raw_expression.tsv \
				./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-T,kUmNA-T" \
		--output ./test6/kUmX-T.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kNmA-T,kNmNA-T" \
		--output ./test6/kNmX-T.tsv


test6_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
				./test6/kUmX-T.tsv \
				./test6/kNmX-T.tsv \
				config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test6/kUmX-T.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test6/fuzzy_kUmX-T/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test6/kNmX-T.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test6/fuzzy_kNmX-T/


test6_report:	generate_report.py \
				./test6/kUmX-T.tsv \
				./test6/kNmX-T.tsv \
				./test6/fuzzy_kUmX-T/concepts_detailed.json \
				./test6/fuzzy_kNmX-T/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test6/ \
		--label1 kUmX-T \
		--label2 kNmX-T \
		--output ./test6/comparison_reports/


test7_dataset_prepare:	prepare_submatrix.py \
						./data/raw_expression.tsv \
						./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-N,kUmNA-N,kNmA-N,kNmNA-N" \
		--output ./test7/data_set/neutrophil.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-T,kUmNA-T,kNmA-T,kNmNA-T" \
		--output ./test7/data_set/T-cell.tsv


test7_dataset_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						./test7/data_set/neutrophil.tsv \
						./test7/data_set/T-cell.tsv \
						config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test7/data_set/neutrophil.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test7/data_set/fuzzy_neutrophil/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test7/data_set/T-cell.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test7/data_set/fuzzy_T-cell/


test7_dataset_report:	generate_report.py \
						./test7/data_set/neutrophil.tsv \
						./test7/data_set/T-cell.tsv \
						./test7/data_set/fuzzy_neutrophil/concepts_detailed.json \
						./test7/data_set/fuzzy_T-cell/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test7/data_set/ \
		--label1 neutrophil \
		--label2 T-cell \
		--output ./test7/data_set/comparison_reports/


test7_normalweight_prepare:	prepare_submatrix.py \
							./data/raw_expression.tsv \
							./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kNmA-N,kNmNA-N" \
		--output ./test7/normalweight/neutrophil.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kNmA-T,kNmNA-T" \
		--output ./test7/normalweight/T-cell.tsv


test7_normalweight_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
							./test7/normalweight/neutrophil.tsv \
							./test7/normalweight/T-cell.tsv \
							config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test7/normalweight/neutrophil.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test7/normalweight/fuzzy_neutrophil/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test7/normalweight/T-cell.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test7/normalweight/fuzzy_T-cell/


test7_normalweight_report:	generate_report.py \
							./test7/normalweight/neutrophil.tsv \
							./test7/normalweight/T-cell.tsv \
							./test7/normalweight/fuzzy_neutrophil/concepts_detailed.json \
							./test7/normalweight/fuzzy_T-cell/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test7/normalweight/ \
		--label1 neutrophil \
		--label2 T-cell \
		--output ./test7/normalweight/comparison_reports/


test7_overweight_prepare:	prepare_submatrix.py \
							./data/raw_expression.tsv \
							./data/metadata.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-N,kUmNA-N" \
		--output ./test7/overweight/neutrophil.tsv
	$(PYTHON) prepare_submatrix.py --data ./data/ \
		--category "kUmA-T,kUmNA-T" \
		--output ./test7/overweight/T-cell.tsv


test7_overweight_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						./test7/overweight/neutrophil.tsv \
						./test7/overweight/T-cell.tsv \
						config.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test7/overweight/neutrophil.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test7/overweight/fuzzy_neutrophil/
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx ./test7/overweight/T-cell.tsv \
		--concept constraints.json \
		--config config.json \
		--threads 8 \
		--output ./test7/overweight/fuzzy_T-cell/


test7_overweight_report:	generate_report.py \
							./test7/overweight/neutrophil.tsv \
							./test7/overweight/T-cell.tsv \
							./test7/overweight/fuzzy_neutrophil/concepts_detailed.json \
							./test7/overweight/fuzzy_T-cell/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test7/overweight/ \
		--label1 neutrophil \
		--label2 T-cell \
		--output ./test7/overweight/comparison_reports/


test7_joint_report:	generate_report.py \
					./test7/overweight/neutrophil.tsv \
					./test7/overweight/T-cell.tsv \
					./test7/normalweight/neutrophil.tsv \
					./test7/normalweight/T-cell.tsv \
					./test7/overweight/fuzzy_neutrophil/concepts_detailed.json \
					./test7/overweight/fuzzy_T-cell/concepts_detailed.json \
					./test7/normalweight/fuzzy_neutrophil/concepts_detailed.json \
					./test7/normalweight/fuzzy_T-cell/concepts_detailed.json
	$(PYTHON) generate_report.py --dir ./test7/ \
		--label1 neutrophil \
		--label2 T-cell \
		--group1 overweight \
		--group2 normalweight \
		--output ./test7/joint_reports/


