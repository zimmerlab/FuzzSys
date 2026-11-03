R = Rscript
PYTHON = python3
FUZZSYS = PATH_TO_FUZZSYS_DIRECTORY
BASE = PATH_TO_BASE_DIRECTORY
# switch here for different data sets
DIR = $(BASE)/miRNA
#DIR = $(BASE)/RNA

DESeq2_log2FC: DESeq2_log2FC_concept DESeq2_log2FC_fuzzy
DESeq2_padj: DESeq2_padj_concept DESeq2_padj_fuzzy
DESeq2: DESeq2_test DESeq2_log2FC DESeq2_padj DESeq2_specific
raw_expression: raw_expression_concept raw_numerator_fuzzy raw_denominator_fuzzy
norm_expression: norm_expression_concept norm_numerator_fuzzy norm_denominator_fuzzy
raw_fuzzy: raw_expression raw_log2FC_concept raw_log2FC_fuzzy raw_log2FC_combine
norm_fuzzy: norm_expression norm_log2FC_concept norm_log2FC_fuzzy norm_log2FC_combine
raw_log2FC: raw_fuzzy raw_expression_flows raw_log2FC_flows raw_log2FC_specific
norm_log2FC: norm_fuzzy norm_expression_flows norm_log2FC_flows norm_log2FC_specific
all: DESeq2 raw_log2FC norm_log2FC
.PHONY: download DESeq2 raw_log2FC norm_log2FC


download:	main_download.R
	$(R) main_download.R $(BASE)
	cp -R ./miRNA/ ./RNA/ $(BASE)


DESeq2_test:	main_DESeq2_test.R \
				$(DIR)/data/raw_counts.tsv \
				$(DIR)/data/metadata.tsv
	$(R) main_DESeq2_test.R $(DIR)/data/raw_counts.tsv $(DIR)/data/metadata.tsv $(DIR)/data/


DESeq2_log2FC_concept:	$(FUZZSYS)/main_concepts.py \
						$(DIR)/data/DESeq2_log2FC.tsv \
						$(DIR)/config/DESeq2_log2FC.json
	$(PYTHON) $(FUZZSYS)/main_concepts.py --mtx $(DIR)/data/DESeq2_log2FC.tsv \
		--config $(DIR)/config/DESeq2_log2FC.json \
		--threads 8 \
		--output $(DIR)/FV_DESeq2/concepts_log2FC/


DESeq2_log2FC_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						$(DIR)/data/DESeq2_log2FC.tsv \
						$(DIR)/FV_DESeq2/concepts_log2FC/concepts_detailed.json \
						$(DIR)/config/DESeq2_log2FC.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/DESeq2_log2FC.tsv \
		--concept $(DIR)/FV_DESeq2/concepts_log2FC/concepts_detailed.json \
		--config $(DIR)/config/DESeq2_log2FC.json \
		--threads 4 \
		--output $(DIR)/FV_DESeq2/log2FC/


DESeq2_padj_concept:	$(FUZZSYS)/main_concepts.py \
						$(DIR)/data/DESeq2_padj.tsv \
						$(DIR)/config/DESeq2_padj.json
	$(PYTHON) $(FUZZSYS)/main_concepts.py --mtx $(DIR)/data/DESeq2_padj.tsv \
		--config $(DIR)/config/DESeq2_padj.json \
		--threads 8 \
		--output $(DIR)/FV_DESeq2/concepts_padj/


DESeq2_padj_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
					$(DIR)/data/DESeq2_padj.tsv \
					$(DIR)/FV_DESeq2/concepts_padj/concepts_detailed.json \
					$(DIR)/config/DESeq2_padj.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/DESeq2_padj.tsv \
		--concept $(DIR)/FV_DESeq2/concepts_padj/concepts_detailed.json \
		--config $(DIR)/config/DESeq2_padj.json \
		--threads 4 \
		--output $(DIR)/FV_DESeq2/padj/


DESeq2_specific:	main_select_DESeq2.py \
					$(DIR)/data/metadata.tsv \
					$(DIR)/data/DESeq2_log2FC.tsv \
					$(DIR)/data/DESeq2_padj.tsv \
					$(wildcard $(DIR)/FV_DESeq2/log2FC/fuzzy_values/fuzzyValues_*.tsv) \
					$(wildcard $(DIR)/FV_DESeq2/padj/fuzzy_values/fuzzyValues_*.tsv)
	$(PYTHON) main_select_DESeq2.py --data $(DIR)/data/ \
		--fuzzy $(DIR)/FV_DESeq2/ \
		--clusterCol context \
		--output $(DIR)/results/


raw_prepare_matrix:	main_preparation.py \
					$(DIR)/data/raw_counts.tsv \
					$(DIR)/data/metadata.tsv \
					$(DIR)/config/preparation.json
	$(PYTHON) main_preparation.py --mtx $(DIR)/data/raw_counts.tsv \
		--metadata $(DIR)/data/metadata.tsv \
		--config $(DIR)/config/preparation.json \
		--output $(DIR)/data/raw_counts/


raw_expression_concept:	$(FUZZSYS)/main_concepts.py \
						$(wildcard $(DIR)/data/raw_counts/denominator_submatrices/submatrix_*.tsv) \
						$(DIR)/config/expression.json
	for f in `ls $(DIR)/data/raw_counts/denominator_submatrices/ | grep -E submatrix_.+\.tsv`; do \
		g="`echo $$f | sed -e 's/\.tsv//g'`"; \
		$(PYTHON) $(FUZZSYS)/main_concepts.py --mtx $(DIR)/data/raw_counts/denominator_submatrices/$$f \
			--config $(DIR)/config/expression.json \
			--threads 12 \
			--output $(DIR)/FV_raw_expression/concepts/$$g/; \
	done


raw_numerator_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						$(wildcard $(DIR)/data/raw_counts/numerator_submatrices/submatrix_*.tsv) \
						$(wildcard $(DIR)/FV_raw_expression/concepts/submatrix_*/concepts_detailed.json) \
						$(DIR)/config/expression.json
	for f in `ls $(DIR)/data/raw_counts/numerator_submatrices/ | grep -E submatrix_.+\.tsv`; do \
		g="`echo $$f | sed -e 's/\.tsv//g'`"; \
		$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/raw_counts/numerator_submatrices/$$f \
			--concept $(DIR)/FV_raw_expression/concepts/$$g/concepts_detailed.json \
			--config $(DIR)/config/expression.json \
			--threads 4 \
			--output $(DIR)/FV_raw_expression/numerator/$$g/; \
	done


raw_denominator_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						$(wildcard $(DIR)/data/raw_counts/denominator_submatrices/submatrix_*.tsv) \
						$(wildcard $(DIR)/FV_raw_expression/concepts/submatrix_*/concepts_detailed.json) \
						$(DIR)/config/expression.json
	for f in `ls $(DIR)/data/raw_counts/denominator_submatrices/ | grep -E submatrix_.+\.tsv`; do \
		g="`echo $$f | sed -e 's/\.tsv//g'`"; \
		$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/raw_counts/denominator_submatrices/$$f \
			--concept $(DIR)/FV_raw_expression/concepts/$$g/concepts_detailed.json \
			--config $(DIR)/config/expression.json \
			--threads 4 \
			--output $(DIR)/FV_raw_expression/denominator/$$g/; \
	done


raw_expression_flows:	main_input_expression.py \
						$(wildcard $(DIR)/FV_raw_expression/numerator/submatrix_*/fuzzy_values/fuzzyValues_*.tsv) \
						$(wildcard $(DIR)/FV_raw_expression/numerator/submatrix_*/evaluations/summary.tsv) \
						$(wildcard $(DIR)/FV_raw_expression/denominator/submatrix_*/fuzzy_values/fuzzyValues_*.tsv) \
						$(wildcard $(DIR)/FV_raw_expression/denominator/submatrix_*/evaluations/summary.tsv) \
						$(DIR)/data/metadata.tsv \
						$(DIR)/config/FlowSets.json
	$(PYTHON) main_input_expression.py --expression $(DIR)/FV_raw_expression/ \
		--metadata $(DIR)/data/metadata.tsv \
		--config $(DIR)/config/FlowSets.json \
		--output $(DIR)/FlowSets_input/raw_counts/


raw_log2FC_concept:	$(FUZZSYS)/main_concepts.py \
					$(DIR)/data/raw_counts/paired_log2FC.tsv \
					$(DIR)/config/data_log2FC.json
	$(PYTHON) $(FUZZSYS)/main_concepts.py --mtx $(DIR)/data/raw_counts/paired_log2FC.tsv \
		--config $(DIR)/config/data_log2FC.json \
		--threads 12 \
		--output $(DIR)/FV_raw_log2FC/


raw_log2FC_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
					$(DIR)/data/raw_counts/paired_log2FC.tsv \
					$(DIR)/FV_raw_log2FC/concepts_detailed.json \
					$(DIR)/config/data_log2FC.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/raw_counts/paired_log2FC.tsv \
		--concept $(DIR)/FV_raw_log2FC/concepts_detailed.json \
		--config $(DIR)/config/data_log2FC.json \
		--threads 4 \
		--output $(DIR)/FV_raw_log2FC/


raw_log2FC_combine:	main_combineFV.py \
					$(wildcard $(DIR)/FV_raw_log2FC/fuzzy_values/*.tsv) \
					$(wildcard $(DIR)/FV_raw_expression/numerator/submatrix_*/fuzzy_values/*.tsv) \
					$(wildcard $(DIR)/FV_raw_expression/denominator/submatrix_*/fuzzy_values/*.tsv) \
					$(DIR)/config/fuzzy_rules.csv
	$(PYTHON) main_combineFV.py --paired $(DIR)/FV_raw_log2FC/fuzzy_values/ \
		--expression $(DIR)/FV_raw_expression/ \
		--rules $(DIR)/config/fuzzy_rules.csv \
		--fuzzySets "NA,-INF,+INF,--,-,o,+,++" \
		--output $(DIR)/FV_raw_log2FC/final_fuzzy_values/


raw_log2FC_flows:	main_input_log2FC.py \
					$(wildcard $(DIR)/FV_raw_log2FC/final_fuzzy_values/fuzzyValues_*.tsv) \
					$(DIR)/data/metadata.tsv \
					$(DIR)/config/FlowSets.json
	$(PYTHON) main_input_log2FC.py --log2FC $(DIR)/FV_raw_log2FC/final_fuzzy_values/ \
		--metadata $(DIR)/data/metadata.tsv \
		--config $(DIR)/config/FlowSets.json \
		--output $(DIR)/FlowSets_input/raw_counts/log2FC.tsv


raw_log2FC_specific:	main_FlowSets.py \
						$(wildcard $(DIR)/FlowSets_input/raw_counts/*.tsv) \
						$(DIR)/config/FlowSets.json
	$(PYTHON) main_FlowSets.py --input $(DIR)/FlowSets_input/raw_counts/ \
		--config $(DIR)/config/FlowSets.json \
		--output $(DIR)/results/raw_counts/


norm_prepare_matrix:	main_preparation.py \
						$(DIR)/data/normalized_counts.tsv \
						$(DIR)/data/metadata.tsv \
						$(DIR)/config/preparation.json
	$(PYTHON) main_preparation.py --mtx $(DIR)/data/normalized_counts.tsv \
		--metadata $(DIR)/data/metadata.tsv \
		--config $(DIR)/config/preparation.json \
		--output $(DIR)/data/norm_counts/


norm_expression_concept:	$(FUZZSYS)/main_concepts.py \
							$(wildcard $(DIR)/data/norm_counts/denominator_submatrices/submatrix_*.tsv) \
							$(DIR)/config/expression.json
	for f in `ls $(DIR)/data/norm_counts/denominator_submatrices/ | grep -E submatrix_.+\.tsv`; do \
		g="`echo $$f | sed -e 's/\.tsv//g'`"; \
		$(PYTHON) $(FUZZSYS)/main_concepts.py --mtx $(DIR)/data/norm_counts/denominator_submatrices/$$f \
			--config $(DIR)/config/expression.json \
			--threads 12 \
			--output $(DIR)/FV_norm_expression/concepts/$$g/; \
	done


norm_numerator_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						$(wildcard $(DIR)/data/norm_counts/numerator_submatrices/submatrix_*.tsv) \
						$(wildcard $(DIR)/FV_norm_expression/concepts/submatrix_*/concepts_detailed.json) \
						$(DIR)/config/expression.json
	for f in `ls $(DIR)/data/norm_counts/numerator_submatrices/ | grep -E submatrix_.+\.tsv`; do \
		g="`echo $$f | sed -e 's/\.tsv//g'`"; \
		$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/norm_counts/numerator_submatrices/$$f \
			--concept $(DIR)/FV_norm_expression/concepts/$$g/concepts_detailed.json \
			--config $(DIR)/config/expression.json \
			--threads 4 \
			--output $(DIR)/FV_norm_expression/numerator/$$g/; \
	done


norm_denominator_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
						$(wildcard $(DIR)/data/norm_counts/denominator_submatrices/submatrix_*.tsv) \
						$(wildcard $(DIR)/FV_norm_expression/concepts/submatrix_*/concepts_detailed.json) \
						$(DIR)/config/expression.json
	for f in `ls $(DIR)/data/norm_counts/denominator_submatrices/ | grep -E submatrix_.+\.tsv`; do \
		g="`echo $$f | sed -e 's/\.tsv//g'`"; \
		$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/norm_counts/denominator_submatrices/$$f \
			--concept $(DIR)/FV_norm_expression/concepts/$$g/concepts_detailed.json \
			--config $(DIR)/config/expression.json \
			--threads 4 \
			--output $(DIR)/FV_norm_expression/denominator/$$g/; \
	done


norm_expression_flows:	main_input_expression.py \
						$(wildcard $(DIR)/FV_norm_expression/numerator/submatrix_*/fuzzy_values/fuzzyValues_*.tsv) \
						$(wildcard $(DIR)/FV_norm_expression/numerator/submatrix_*/evaluations/summary.tsv) \
						$(wildcard $(DIR)/FV_norm_expression/denominator/submatrix_*/fuzzy_values/fuzzyValues_*.tsv) \
						$(wildcard $(DIR)/FV_norm_expression/denominator/submatrix_*/evaluations/summary.tsv) \
						$(DIR)/data/metadata.tsv \
						$(DIR)/config/FlowSets.json
	$(PYTHON) main_input_expression.py --expression $(DIR)/FV_norm_expression/ \
		--metadata $(DIR)/data/metadata.tsv \
		--config $(DIR)/config/FlowSets.json \
		--output $(DIR)/FlowSets_input/norm_counts/


norm_log2FC_concept:	$(FUZZSYS)/main_concepts.py \
						$(DIR)/data/norm_counts/paired_log2FC.tsv \
						$(DIR)/config/data_log2FC.json
	$(PYTHON) $(FUZZSYS)/main_concepts.py --mtx $(DIR)/data/norm_counts/paired_log2FC.tsv \
		--config $(DIR)/config/data_log2FC.json \
		--threads 12 \
		--output $(DIR)/FV_norm_log2FC/


norm_log2FC_fuzzy:	$(FUZZSYS)/main_fuzzify.py \
					$(DIR)/data/norm_counts/paired_log2FC.tsv \
					$(DIR)/FV_norm_log2FC/concepts_detailed.json \
					$(DIR)/config/data_log2FC.json
	$(PYTHON) $(FUZZSYS)/main_fuzzify.py --mtx $(DIR)/data/norm_counts/paired_log2FC.tsv \
		--concept $(DIR)/FV_norm_log2FC/concepts_detailed.json \
		--config $(DIR)/config/data_log2FC.json \
		--threads 4 \
		--output $(DIR)/FV_norm_log2FC/


norm_log2FC_combine:	main_combineFV.py \
						$(wildcard $(DIR)/FV_norm_log2FC/fuzzy_values/*.tsv) \
						$(wildcard $(DIR)/FV_norm_expression/numerator/submatrix_*/fuzzy_values/*.tsv) \
						$(wildcard $(DIR)/FV_norm_expression/denominator/submatrix_*/fuzzy_values/*.tsv) \
						$(DIR)/config/fuzzy_rules.csv
	$(PYTHON) main_combineFV.py --paired $(DIR)/FV_norm_log2FC/fuzzy_values/ \
		--expression $(DIR)/FV_norm_expression/ \
		--rules $(DIR)/config/fuzzy_rules.csv \
		--fuzzySets "NA,-INF,+INF,--,-,o,+,++" \
		--output $(DIR)/FV_norm_log2FC/final_fuzzy_values/


norm_log2FC_flows:	main_input_log2FC.py \
					$(wildcard $(DIR)/FV_norm_log2FC/final_fuzzy_values/fuzzyValues_*.tsv) \
					$(DIR)/data/metadata.tsv \
					$(DIR)/config/FlowSets.json
	$(PYTHON) main_input_log2FC.py --log2FC $(DIR)/FV_norm_log2FC/final_fuzzy_values/ \
		--metadata $(DIR)/data/metadata.tsv \
		--config $(DIR)/config/FlowSets.json \
		--output $(DIR)/FlowSets_input/norm_counts/log2FC.tsv


norm_log2FC_specific:	main_FlowSets.py \
						$(wildcard $(DIR)/FlowSets_input/norm_counts/*.tsv) \
						$(DIR)/config/FlowSets.json
	$(PYTHON) main_FlowSets.py --input $(DIR)/FlowSets_input/norm_counts/ \
		--config $(DIR)/config/FlowSets.json \
		--output $(DIR)/results/norm_counts/


