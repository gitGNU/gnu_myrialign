
#
#    Copyright 2008 Paul Harrison
#
#    This file is part of Myrialign.
#    
#    Myrialign is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    Myrialign is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
#

"""

    Analyse the output of "myr align" in various ways.

"""

import sys, numpy, os.path

import sequence

class Error(Exception): pass
class Bad_option(Error): pass

class Hit:
    """ Structure for storing a hit. Members:
    
        name
	direction (read direction 'fwd'/'rev')
        start 
	end
	read_ali
	ref_ali """

def get_option(argv, option):
    argv = argv[:]
    has_option = False
    while True:
        try:
	    location = argv.index(option)
	except ValueError: #Not found
	    break
        has_option = True
	del argv[location]

    return has_option, argv

def get_option_value(argv, option, conversion_function, default):
    argv = argv[:]
    value = default
    while True:
        try:
	    location = argv.index(option)
	except ValueError: #Not found
	    break
	    
	if location == len(argv)-1 :
	    raise Bad_option('Option %s requires a paramter' % option)
	
	try:
	    value = conversion_function(argv[location+1])
	except Exception:
	    raise Bad_option('Option for %s not in expected format' % option)
	
	del argv[location:location+2]

    return value, argv


def show_default_options():
     print >> sys.stderr, '    -u    - Only count reads that have a unique hit'
     print >> sys.stderr, ''
     print >> sys.stderr, 'Alignments are of lower quality at the start and end. Clipping' 
     print >> sys.stderr, 'these will produce higher quality results.'
     print >> sys.stderr, ''
     print >> sys.stderr, '    -s n  - Clip n bases from read start'
     print >> sys.stderr, '    -e n  - Clip n bases from read end'

def clip_alignment(read_ali, ref_ali, clip_start, clip_end):
    n_start = 0
    while read_ali and clip_start > 0:
        if read_ali[0] != '-':
	    clip_start -= 1
	read_ali = read_ali[1:]
	ref_ali = ref_ali[1:]
	n_start += 1

    n_end = 0
    while read_ali and clip_end > 0:
        if read_ali[-1] != '-':
	    clip_end -= 1
	read_ali = read_ali[:-1]
	ref_ali = ref_ali[:-1]
	n_end += 1
    
    return read_ali, ref_ali, n_start, n_end


class Table:
    def __init__(self, members):
	self.store_size = 1
	self.length = 0
	self.members = [ ]
	self.indicies = { }
        for name, dtype in members:
	    self.members.append(name)
	    setattr(self, name, numpy.empty(self.store_size, dtype))
	
    def resize(self, length):
	self.length = length
	self.indicies = { }
	if length > self.store_size:	
	    self.store_size = length*5//4
	    for name in self.members:
	        getattr(self, name).resize(self.store_size)

    def iter_groups(self, name):
        column = getattr(self, name)
        order = numpy.argsort(column)
	start = 0
	while start < len(order):
	    end = start + 1
	    while end < len(order) and column[order[start]] == column[order[end]]:
	        end += 1
	    yield order[start:end]
	    start = end
	    

def read_files(argv):
    clip_start, argv = get_option_value(argv, '-s', int, 0)
    clip_end, argv = get_option_value(argv, '-e', int, 0)

    if len(argv) < 2:
        raise Bad_option('Expected at least two filenames, a reference genome and output from myr align')

    reference = sequence.sequence_file_iterator(argv[0]).next()[1] 

    #read_hits = { }
    
    hits = Table((
        ('name', 'object'),
	('forward', 'bool'),
	('read_ali', 'object'),
	('ref_ali', 'object'),
	('start', 'int32'),
	('end', 'int32'),
	('n_errors', 'int32'),
    ))

    nth = 0
    for filename in argv[1:]:
	for line in open(filename,'rb'):
            if not line.endswith('\n'): continue
	    if line.startswith('#'): continue

	    #hit = Hit()
	    #hit.name, hit.direction, hit.n_errors, span, hit.read_ali, hit.ref_ali = line.rstrip().split()
	    i = hits.length
	    hits.resize(i+1)
	    hits.name[i], direction, n_errors, span, hits.read_ali[i], hits.ref_ali[i] = line.rstrip().split()
	    start, end = span.split('..')
	    hits.start[i] = int(start)-1
	    hits.end[i] = int(end)
	    hits.n_errors[i] = int(n_errors)	    
	    hits.forward[i] = (direction == 'fwd')

            if clip_start or clip_end:
	        if hits.forward[i]:
		    hits.read_ali[i], hits.ref_ali[i], clipped_start, clipped_end = clip_alignment(hits.read_ali[i], hits.ref_ali[i], clip_start, clip_end)
	        else:
		    hits.read_ali[i], hits.ref_ali[i], clipped_start, clipped_end = clip_alignment(hits.read_ali[i], hits.ref_ali[i], clip_end, clip_start)
	        hits.start[i] += clipped_start
	        hits.end[i] -= clipped_end

	    #if hits.name[i] not in read_hits: 
        	#read_hits[hit.name] = [ ]
	    #read_hits[hit.name].append(hit) 

	    nth += 1
	    if nth % 10000 == 0:
        	sys.stderr.write('Loading hits: %d            \r' % nth)
		sys.stderr.flush()
    
    return reference, hits




def artplot(argv):
    try:
        only_single, argv = get_option(argv, '-u')
        reference, hits = read_files(argv)
    except Bad_option, error:
        print >> sys.stderr, ''
	print >> sys.stderr, 'myr artplot [options] <reference genome> <myr align output> [<myr align output>...]'
	print >> sys.stderr, ''
	print >> sys.stderr, 'Options:'
	print >> sys.stderr, ''
	show_default_options()
	print >> sys.stderr, ''
	print >> sys.stderr, error[0]
	return 1

    size = len(reference)

    coverage = numpy.zeros(size, 'float64')
    insertions = numpy.zeros(size, 'float64')
    deletions = numpy.zeros(size, 'float64')
    substitutions = numpy.zeros(size, 'float64')
    base_counts = numpy.zeros((size,5), 'float64')
    base_map = {'A':0,'T':1,'C':2,'G':3,'N':4}
    
    nth = 0
    for group in hits.iter_groups('name'):
	if only_single and len(group) > 1: continue

	weight = 1.0 / len(group)

	for i in group:    
	    pos = hits.start[i]
	    read_ali = hits.read_ali[i]
	    ref_ali = hits.ref_ali[i]
	    for j in xrange(len(read_ali)):
		a = read_ali[j]
		b = ref_ali[j]

        	if a == '-':
		    deletions[pos] += weight
		    coverage[pos] += weight
		    pos += 1
		elif b == '-':
		    insertions[pos] += weight
		else:
		    if a != b:
	                substitutions[pos] += weight

	            coverage[pos] += weight
		    base_counts[pos, base_map[a]] += weight
		    pos += 1

	    nth += 1
	    if nth % 10000 == 0:
        	sys.stderr.write('Processing: %d          \r' % nth)
		sys.stderr.flush()

    sys.stderr.write(' %d hits\n' % nth)

    def save(filename, array):
	print 'Writing', filename
	f = open(filename, 'wb')
	for i in array:
            print >> f, i
	f.close()

    prefix = os.path.splitext(os.path.basename(argv[0]))[0] + '-'

    save(prefix+'coverage.txt', coverage)
    save(prefix+'insertions.txt', insertions)
    save(prefix+'deletions.txt', deletions)
    save(prefix+'substitutions.txt', substitutions)
    
    base_counts = base_counts[:,:4]
    base_total = numpy.sum(base_counts, 1)
    surprise = numpy.zeros((size,4),'float64')
    for i in xrange(4):
        good = base_counts[:,i] > 0
        surprise[good,i] = numpy.log( base_counts[good,i] / base_total[good] ) / numpy.log(0.5)
    entropy = numpy.zeros(size, 'float64')
    good = base_total > 0
    entropy[good] = numpy.sum(base_counts[good,:] * surprise[good,:], 1) / base_total[good]
    save(prefix+'confusion.txt', entropy)


def textdump(argv):
    try:
        only_single, argv = get_option(argv, '-u')
        reference, hits = read_files(argv)
    except Bad_option, error:
        print >> sys.stderr, ''
	print >> sys.stderr, 'myr textdump [options] <reference genome> <myr align output> [<myr align output>...]'
	print >> sys.stderr, ''
	print >> sys.stderr, 'Options:'
	print >> sys.stderr, ''
	show_default_options()
	print >> sys.stderr, ''
	print >> sys.stderr, error[0]
	return 1

    size = len(reference)

    todo = { }
    for group in hits.iter_groups('name'):
        if only_single and len(group) > 1: continue	
	
	for i in group:
	    if hits.start[i] not in todo: 
	        todo[hits.start[i]] = []
	    todo[hits.start[i]].append(i)

    lanes = [ ]
    def find_lane():
        for start in (0,4,2,6,1,5,3,7):
            for i in xrange(start,len(lanes),8):
	        if lanes[i] is None: 
	            return i
	lanes.extend([None]*8)
	return find_lane()
    
    pad = ' '*5

    total_with_a_hit = 0    
    for pos, ref_nuc in enumerate(sequence.string_from_sequence(reference)):    
        while lanes and lanes[-1] is None:
	    del lanes[-1]
    
        if pos in todo:
	    for i in todo[pos]:
	        lane_no = find_lane()
		lanes[lane_no] = [hits.ref_ali[i]+pad,hits.read_ali[i]+pad]

        to_show = [ ref_nuc ]
	
	for i in xrange(len(lanes)):
	    if lanes[i] is None:
	        to_show.append('')
		continue
	
	    lane = lanes[i]
            n = 0
	    while n < len(lane[0])-1 and lane[0][n] == '-': n += 1
	    n += 1
	    to_show.append(lane[1][:n])
	    
	    if n == len(lane[0]):
	        lanes[i] = None		
	    else:
	        lane[0] = lane[0][n:]
	        lane[1] = lane[1][n:]

        counts = { }
	for item in to_show[1:]:
	    if item and item != ' ':
	        counts[item] = counts.get(item,0)+1

        consensus = None
	for item in counts:
	    if not consensus or counts[item] > counts[consensus]: #TODO: what if equal?
	        consensus = item
		
        interesting = False
	confusing = False
	for item in counts:
	    if counts[item] > 5: #Arbitrary...
	        if item != to_show[0]:
		    interesting = True
		if item != consensus:
		    confusing = True
	
	if counts:
	    total_with_a_hit += 1
		
        if 0:
            counts = { }
	    for item in to_show[1:]:
		if item and item != ' ':
	            counts[item] = counts.get(item,0)+1
	    unique_items = counts.keys()
	    unique_items.sort(lambda a,b: cmp(counts[b],counts[a]))

	    print to_show[0],
	    for item in unique_items:
		print item + 'x' + str(counts[item]),
            print
        else:
            max_to_show = max([ len(item) for item in to_show ])
	    max_to_show = max(2, max_to_show)
	    to_show = [ ' '*(max_to_show-len(item)) + item for item in to_show ]
	    #print to_show
	    for i in xrange(max_to_show):
		row = [ item[i] for item in to_show ]
		if i == max_to_show-1:
		    print '%9d'%(pos+1),
		else:
		    print ' '*9,
		print row[0], ((confusing and ' ?') or (interesting and '! ') or '  '), ''.join(row[1:])

    print >> sys.stderr, 'Proportion of reference with at least one hit: %.2f%%' % ( 100.0*float(total_with_a_hit)/len(reference) )
