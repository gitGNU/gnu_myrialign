
import random, os, sys

import cache, sequence, align

#TODO: should look at datestamp on filenames

def sample(working_dir, read_files, n_samples):
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
    
sample = cache.cached(sample)


def invoke_align(working_dir, reference_filename, read_filename):
    print >> sys.stderr, 'Aligning'
    #Hmm
    old_stdout = sys.stdout
    sys.stdout = open(os.path.join(working_dir,'hits.myr'), 'wb')

    try:
        assert align.main(['5','1',reference_filename,read_filename]) == 0
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

invoke_align = cache.cached(invoke_align)


def main(argv):
    if len(argv) < 2:
        print >> sys.stderr, ''
	print >> sys.stderr, 'myr assess <contigs file> <reads> [<reads> ...]'
	print >> sys.stderr, ''
	return 1
    
    sample_dir = sample(argv[1:], 100)    
    sample_file = os.path.join(sample_dir, 'sample.fna')

    hit_dir = invoke_align(argv[0], sample_file)
    hit_file = os.path.join(hit_dir, 'hits.myr')

    hits = [ ]
    for item in sequence.sequence_file_iterator(sample_file):
        hits[item[0]] = [ ]

    for line in open(hit_file, 'rb'):
        line = line.strip()
	if line.startswith('#'): continue

	name, direction, n_errors, span, read_ali, ref_ali = line.rstrip().split()
        hits[name].append(int(n_errors))
    
    n_ambiguous = 0
    n_unhit = 0
    error_count = [ 0 ] * 6 #TODO: make max errors an option
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

