
import random, os, sys

import cache, sequence, align

def sample(read_files, n_samples):
    read_filesigs = [ cache.file_signature(filename) for filename in read_files ]
    read_files = [ item[0] for item in read_filesigs ]

    def callback(working_dir):    
        print >> sys.stderr, 'Sampling'    
	samples = [ ]
	n = 0
	for item in sequence.sequence_files_iterator(read_files):
            n += 1
            if len(samples) < n_samples:
		samples.append(item)
	    elif random.random()*n_samples < n:
		samples[random.randrange(n_samples)] = item

	outfile = open(os.path.join(working_dir,'sample.fna'),'wb')
	for item in samples:
            print >> outfile, '>%s' % item[0]
	    print >> outfile, '%s' % sequence.string_from_sequence(item[1])
	    
    result_dir = cache.get(('assess','sample',n_samples,read_filesigs), callback)
    return os.path.join(result_dir, 'sample.fna')

def invoke_align(reference_filename, read_filename, max_errors):
    reference_filesig = cache.file_signature(reference_filename)
    reference_filename = reference_filesig[0]
    read_filesig = cache.file_signature(read_filename)
    read_filename = read_filesig[0]

    def callback(working_dir):
	print >> sys.stderr, 'Aligning'
	#Hmm
	old_stdout = sys.stdout
	sys.stdout = open(os.path.join(working_dir,'hits.myr'), 'wb')

	try:
	    assert align.main([str(max_errors),'1',reference_filename,read_filename]) == 0
	finally:
	    sys.stdout.close()
	    sys.stdout = old_stdout

    return os.path.join(
        cache.get(('assess','invoke_align1',reference_filesig,read_filesig,max_errors),callback), 
        'hits.myr')

def main(argv):
    if len(argv) < 2:
        print >> sys.stderr, ''
	print >> sys.stderr, 'myr assess <contigs file> <reads> [<reads> ...]'
	print >> sys.stderr, ''
	return 1

    max_errors = 3
    
    sample_file = sample(argv[1:], 1000)    

    hit_file = invoke_align(argv[0], sample_file, max_errors)

    hits = { }
    for item in sequence.sequence_file_iterator(sample_file):
        hits[item[0]] = [ ]

    for line in open(hit_file, 'rb'):
        line = line.strip()
	if line.startswith('#'): continue

	name, direction, n_errors, span, read_ali, ref_ali = line.rstrip().split()
        hits[name].append(int(n_errors))
    
    n_ambiguous = 0
    n_unhit = 0
    error_count = [ 0 ] * (max_errors+1)
    for name in hits:
        if not hits[name]:
	    n_unhit += 1
	elif len(hits[name]) > 1:
	    n_ambiguous += 1
	else:
	    error_count[hits[name][0]] += 1
    
    print 'Sampled', len(hits), 'reads'
    print n_ambiguous, 'hit multiple contigs'
    print n_unhit, 'hit nothing'
    for i in xrange(len(error_count)):
        print '%3d errors: %d' % (i,error_count[i])

