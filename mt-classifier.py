#!/usr/bin/env python

import getopt, sys, re, os, glob
from classifier import tree, NGclassify, consts, datatypes, parse_mhcs
from bioinf.seqs import SeqList

# folder where to find data for haplogroup classification and functional annotation
data_file = os.path.dirname(sys.argv[0])

def usage_old():
	print """\nAssigns haplogroup to contigs and performs functional annotation
		Options:
		-i		Contig file [mtDNAassembly-Contigs.fasta]
		-g		GMAP executable PATH [/usr/local/bin/gmap]
		-D		GMAP mt sequences database location [/usr/local/share/gmapdb]
		-m		GMAP mt sequences database [mt_mhcss]
		-t		GMAP threads [2]
		-b		basename for output files
		"""

def usage():
	print """\nAssigns haplogroup to contigs and performs functional annotation
		Options:
		-i		Contig file [mtDNAassembly-Contigs.fasta]
		-m		MUSCLE executable PATH [/usr/local/bin/muscle]
		-b		basename for output files
		-s		file with most reliable haplogroup prediction
		"""

def parse_gmapf9_line(line):
	parsed = line.split('\t')
	last_field = re.findall(r"[\w']+", parsed[2])
	seq_nuc = parsed[1].partition(' ')[2]
	seq_index = parsed[1].partition(' ')[0]
	ref_pos = int(last_field[1])
	ref_nuc = parsed[2][-1]
	return ref_pos, ref_nuc, seq_nuc, seq_index

def parse_gmapf9_file(inhandle):
	contigs_mappings = [[]]
	h = inhandle.readlines()
	c = 0
	mutations = []
	while c < len(h):
		# end coordinate of last contig
		if c == len(h)-1:
			contigs_mappings[-1].append(parse_gmapf9_line(h[c])[0])
		if h[c][0] != '>':
			ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
			# insertion
			if ref_nuc == ' ' and seq_nuc != ' ':
				# gmap assigns the position of the next nucleotide to the insertion
				pos_ins = ref_pos - 1
				ins = [seq_nuc]
				c += 1
				ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
				while c < len(h) and (ref_nuc == ' ' and seq_nuc != ' '):
					ins.append(seq_nuc)
					c += 1
					ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
				mut = datatypes.Insertion("%d.%s" % (pos_ins, ''.join(ins)))
				mutations.append(mut)
				#print "%d.%s" % (pos_ins, ''.join(ins))
			# deletion
			elif ref_nuc != ' ' and seq_nuc == ' ':
				pos_del = ref_pos
				c += 1
				ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
				while c < len(h) and (ref_nuc != ' ' and seq_nuc == ' '):
					c += 1
					ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
				if pos_del == ref_pos-1:
					print "%dd" % (pos_del)
					mut = datatypes.Deletion("%dd" % pos_del)
					mutations.append(mut)
				else:
					print "%d-%dd" % (pos_del, ref_pos-1)
					mut = datatypes.Deletion("%d-%dd" % (pos_del, ref_pos-1))
					mutations.append(mut)
			# mismatch
			elif ref_nuc != seq_nuc:
				if seq_nuc != 'N':
					# Transition
					if (ref_nuc in consts.PUR and seq_nuc in consts.PUR) or (ref_nuc in consts.PYR and seq_nuc in consts.PYR):
						print "%d%s" % (ref_pos, seq_nuc)
						mut = datatypes.Transition(ref_pos)
						mutations.append(mut)
					# Transversion
					if (ref_nuc in consts.PUR and seq_nuc in consts.PYR) or (ref_nuc in consts.PYR and seq_nuc in consts.PUR):
						mut = datatypes.Transversion("%d%s" % (ref_pos, seq_nuc))
						mutations.append(mut)
				c += 1
			else:
				c += 1
		else:
			# first contig
			if len(contigs_mappings) == 1 and len(contigs_mappings[-1]) == 0:
				contigs_mappings[-1].append(parse_gmapf9_line(h[c+1])[0])
			# all the others
			else:
				contigs_mappings[-1].append(parse_gmapf9_line(h[c-1])[0])
				contigs_mappings.append([parse_gmapf9_line(h[c+1])[0]])
			c += 1
	# don't know if contig coordinate sorting is needed but I'll do anyway
	contigs_mappings.sort()
	return mutations, contigs_mappings

def merge_tables(f, g, h):
    fgh = f + g + h
    mergedlist = []
    for jj in fgh:
        if jj not in mergedlist:
            mergedlist.append(jj)
    o = []
    o.append(["", "RSRS", "MHCS", "rCRS"])
    y = "yes"
    n = ""
    for i in mergedlist:
        if i in f and i in g and i in h:
            o.append([i.pprint(),y,y,y])
        elif i in f and i in g:
            o.append([i.pprint(),y,y,n])
        elif i in f and i in h:
            o.append([i.pprint(),y,n,y])
        elif i in g and i in h:
            o.append([i.pprint(),n,y,y])
        elif i in f:
            o.append([i.pprint(),y,n,n])
        elif i in g:
            o.append([i.pprint(),n,y,n])
        elif i in h:
            o.append([i.pprint(),n,n,y])
    return o

def align_sequence(muscle_exe, sequence, rif=None, ):
    """sequence is a datatypes.Sequence, rif"""
    if rif is None:
        rif = datatypes.Sequence('RSRS', consts.RCRS)
    seq_diff = NGclassify.SequenceDiff()
    #print "Aligning sequence %s" % sequence.name
    seq_diff.gen_diff(muscle_exe, rif, datatypes.Sequence(sequence.name, str(sequence)))
    #print "-"*30
    return seq_diff

def h_analysis(htrees, seq_diff, regions, mhcs_dict):
    a = NGclassify.Classify()
    #print "Classification of sequence %s" % seq_diff.obj.name
    for htree, name in htrees:
        print "Classification according to tree:", name
        a.classify_by_tree(htree, seq_diff, regions)
        #print "start is ", seq_diff.start
        #print "end is ", seq_diff.end
        #print "haplo_stats: ", a.haplo_stats
        print "genome_state is ", a.get_genome_state()
        (haplo_stats_sorted, haplo_best) = a.prediction_sorting()
        print haplo_best
        #print "haplo_stats_sorted is:\n", haplo_stats_sorted
        print "="*20
        mhcss = a.get_mhcss(mhcs_dict)
    print '-'*30
    return a

def load_sequences(fname):
    a = SeqList()
    a.load_file(fname)
    print "Loaded %d contig sequences" % len(a)
    return a

def write_output(class_obj, seq_diff, seq_diff_mhcs, seq_diff_rcrs, merged_tables, outfile):
    print "Writing results for sequence %s" % outfile
    class_obj.pprint(open(outfile + '.csv', 'w'))
    class_obj.pprint_sorted(open(outfile + '.sorted.csv', 'w'))
    merged_tables_file = open(outfile + '_merged_diff.csv', 'w')
    for row in merged_tables:
        merged_tables_file.write(','.join(row)+'\n')

def main_mt_hpred():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:m:b:s:")
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit()
	#print opts, args
	contig_file = 'mtDNAassembly-contigs.fasta'
	muscle_exe='/usr/local/bin/muscle'
	basename='mtDNAassembly-contigs'
	best_results_file = 'mt_classification_best_results.csv'
	#print opts
	for o,a in opts:
		#print "option", o, "argument", a
		if o == "-h":
			usage()
			sys.exit()
		elif o == "-i": contig_file = a
		elif o == "-m": muscle_exe = a
		elif o == "-b": basename = a
		elif o == "-s": best_results_file = a
		else:
			assert False, "Unhandled option."

	print "Your best results file is ", best_results_file
	# sample name
	f = os.path.abspath(contig_file)
	#sample_name = f.split('/')[-2].split('_')[-1]
	sample_name = contig_file.split('-')[0]
	
	# haplogroup tree parsing
	htrees = [(tree.HaplogroupTree(pickle_data=open(data_file + '/phylotree_r15.pickle', 'rb').read()), data_file + '/data/phylotree_r15.pickle')]
	# mhcs parsing
	mhcs_dict = parse_mhcs.parse2mhcs_dict(data_file + '/data/mhcs.tab')
	
	print "\nLoading contig sequences from file %s" % contig_file
	contig_array = load_sequences(contig_file)
	contig_array_seqdiff = [] # lista di liste
	contig_total_seqdiff = [] # lista di varianti
	contig_array_mappings = []
	
	print "\nAligning Contigs to mtDNA reference genome...\n"
	
	# update each contig's SeqDiff
	for x,contig in enumerate(contig_array):
		if x == 0:
			contig_seq_diff = align_sequence(muscle_exe, contig)
			contig_seq_diff.find_segment() # avoid having long gaps at 5' and 3' (not actual gaps but due to the alignment)
			contig_seq_diff.regions.append([contig_seq_diff.start, contig_seq_diff.end])
		else:
			incoming_seqdiff = align_sequence(muscle_exe, contig)
			incoming_seqdiff.find_segment()
			contig_seq_diff.diff_list.extend(incoming_seqdiff.diff_list)
			contig_seq_diff.regions.append([incoming_seqdiff.start, incoming_seqdiff.end])

	print "\nSequence haplogroup assignment\n"
	seq_classify = h_analysis(htrees, contig_seq_diff, contig_seq_diff.regions, mhcs_dict)
	seq_classify.sample_name = sample_name

	#print "\nSequence functional annotation\n"
	print "Contig alignment to MHCS and rCRS"
	m = list(seq_classify.mhcss)[0]
	print "Aligning contigs to MHCS SeqDiff object"
	its_mhcs = datatypes.Sequence(m, mhcs_dict[m])
	#contig_mhcs_total_seqdiff = []
	for x, contig in enumerate(contig_array):
		if x == 0:
			contig_mhcs_seq_diff = align_sequence(muscle_exe, contig, its_mhcs)
			contig_mhcs_seq_diff.find_segment()
			contig_mhcs_seq_diff.regions.append([contig_seq_diff.start, contig_seq_diff.end])
		else:
			incoming_mhcs_seqdiff = align_sequence(muscle_exe, contig, its_mhcs)
			incoming_mhcs_seqdiff.find_segment()
			contig_mhcs_seq_diff.diff_list.extend(incoming_mhcs_seqdiff.diff_list)
			contig_mhcs_seq_diff.regions.append([incoming_mhcs_seqdiff.start, incoming_mhcs_seqdiff.end])

	print "rCRS SeqDiff object"
	rcrs = datatypes.Sequence('rCRS', consts.rcrs)
	for x, contig in enumerate(contig_array):
		if x == 0:
			contig_rcrs_seq_diff = align_sequence(muscle_exe, contig, rcrs)
			contig_rcrs_seq_diff.find_segment()
			contig_rcrs_seq_diff.regions.append([contig_seq_diff.start, contig_seq_diff.end])
		else:
			incoming_rcrs_seqdiff = align_sequence(muscle_exe, contig, rcrs)
			incoming_rcrs_seqdiff.find_segment()
			contig_rcrs_seq_diff.diff_list.extend(incoming_rcrs_seqdiff.diff_list)
			contig_rcrs_seq_diff.regions.append([incoming_rcrs_seqdiff.start, incoming_rcrs_seqdiff.end])

	print "Merging seq_diffs..."
	mergedtables = merge_tables(contig_seq_diff.diff_list, contig_mhcs_seq_diff.diff_list, contig_rcrs_seq_diff.diff_list)
	write_output(seq_classify, contig_seq_diff.diff_list, contig_mhcs_seq_diff.diff_list, contig_rcrs_seq_diff.diff_list, mergedtables, basename)
	open(os.path.join('../', best_results_file), 'a').write(','.join([basename, ';'.join([i[0] for i in seq_classify.haplo_best.items()])])+'\n')
	
if __name__ == "__main__":
	main_mt_hpred()

