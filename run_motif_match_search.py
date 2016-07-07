import motifMatchSearch
import sys
import argparse

if __name__ == "__main__":

	mms = motifMatchSearch.motif_query_handler()
	p_value = '10e-10'

	parser = argparse.ArgumentParser(description='Find motifs in input')
	parser.add_argument('-f', '--fasta')
	parser.add_argument('-b', '--bed')
	parser.add_argument('-bgs', '--background-summary')
	parser.add_argument('-o', '--outpath')
	parser.add_argument("-pfx", '--prefix')
	parser.add_argument('-p', '--p-value')
	parser.add_argument('-fl', '--flank')
	parser.add_argument('-spi', '--one-score-per-input')
	parser.add_argument('-bg', '--print-bg', action='store_true')
	parser.add_argument('-hits', '--print-hits', action='store_true')
	parser.add_argument('-sum', '--print-summary', action='store_true')

	args = parser.parse_args()
	print args

	## set the output type
	out_bg = False
	out_hits = False
	out_summary = False

	if args.print_bg:
		out_bg = True

	if args.print_hits:
		out_hits = True

	if args.print_summary:
		out_summary = True

	if not out_bg and not out_hits and not out_summary:
		print "please select at least one output format (-bg, -hits, -sum)"
		sys.exit(0)

	mms.set_out_bedgraph(out_bg)
	mms.set_out_hits(out_hits)

	## set the output directory
	mms.set_out_path(args.outpath)

	## set the input file
	pfx = "rna_compete_motif"
	mms.set_out_prefix(pfx)

	if args.fasta:
		mms.set_fasta(args.fasta)
		mms.set_out_prefix(args.fasta[0:-4])

	elif args.bed:
		mms.set_bed(args.bed)
		mms.set_out_prefix(args.bed[0:-4])
	else:
		print "Please provide a set of input sequences with -b or -f"
		sys.exit(0)

	## Set the p value threshold for hits
	if args.p_value:
		p_value = args.p_value

	## this controls whether all motifs are counted or a maximum of one per input sequence
	spi = False
	if args.one_score_per_input:
		spi = True
	mms.set_one_count_per_input(spi)

	## if the output is bedgraph, pad it by this amount
	f_range = 75
	mms.set_flank_range(f_range)

	## hardcoded background for now
	background_summary = 'background_annotation_bed/refGene_refGene_exon_plus_flank_unique.summary_bnf.txt'
	mms.set_background_from_summary(background_summary)

	## set the pval threshold for calling hits
	mms.set_p_thresh(float(p_value))

	## run the search and do the summary if requested
	mms.run_search()
	if out_summary:
		mms.summarize()


