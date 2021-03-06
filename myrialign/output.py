
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
#    along with Myrialign.  If not, see <http://www.gnu.org/licenses/>.
#

"""

    Analyse the output of "myr align" in various ways.

"""

import sys, numpy, os.path, sets, heapq

import sequence, sort

class Error(Exception): pass
class Bad_option(Error): pass
class Not_found(Error): pass
class Out_of_bounds(Error): pass

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
    def __init__(self):
	self.store_size = 0
	self.length = 0
	self.indicies = { }
        for name, dtype in self.MEMBERS:
	    setattr(self, name, numpy.empty(self.store_size, dtype))
	
    def resize(self, length):
	self.length = length
	self.indicies = { }
	if length > self.store_size:
	    old_size = self.store_size
	    self.store_size = length*5//4
	    for name, dtype in self.MEMBERS:
	        column = numpy.empty(self.store_size, dtype)
		column[:old_size] = self.__dict__[name]
	        setattr(self,name,column)

    def new_id(self):
        result = self.length
	self.resize(result+1)
	return result

    def you_are_dirty(self):
        if self.indicies:
	    self.indicies = { }

    def index(self, name):
        if name not in self.indicies:
	    self.indicies[name] = numpy.argsort(self.__dict__[name][:self.length])
	return self.indicies[name]

    def _find(self, name, value):
        column = self.__dict__[name]
        order = self.index(name)
	
	lo = 0
	hi = len(order)
	while lo < hi:
            mid = (lo+hi)//2
            if column[order[mid]] < value: lo = mid+1
            else: hi = mid

        return column, order, lo

    def find(self, name, value):
        column, order, i = self._find(name, value)
	if column[order[i]] != value:
	    raise Not_found(self,name,value)
	return order[i]
    
    def find_all(self, name, value):
        column, order, i = self._find(name, value)
	j = i
	while j < len(order) and column[order[j]] == value: 
	    j += 1
	return order[i:j]

    def iter_groups(self, name):
        column = self.__dict__[name]
        order = self.index(name)
	start = 0
	while start < len(order):
	    end = start + 1
	    while end < len(order) and column[order[start]] == column[order[end]]:
	        end += 1
	    yield order[start:end]
	    start = end

class Hits(Table):
    MEMBERS = (
        ('name', 'object'),
	('forward', 'bool'),
	('read_ali', 'object'),
	('ref_ali', 'object'),
	('start', 'int32'),
	('end', 'int32'),
    )   


def iter_hit_file_myrialign(filename):
    ref_name = None
    nth = 0
    for line in open(filename, 'rb'):
        if not line.endswith('\n'): break #Alignment file truncated or still being written
    
        if line.startswith('#'):
            if line.startswith('#Reference:'):
	        ref_name = line.split()[1]
	    continue
	
	name, direction, n_errors, span, read_ali, ref_ali = line.rstrip().split()
	start, end = span.split('..')
	start = int(start)-1
	end = int(end)
	forward = (direction == 'fwd')
	
	nth += 1
	if nth % 1000 == 0:
            sys.stderr.write('Loading hits: %d            \r' % nth)
	    sys.stderr.flush()
	
	yield ref_name, name, forward, start, end, read_ali, ref_ali

def iter_hit_file_maf(filename):
    # BLAT only, for now
    seqs = [ ]
    f = open(filename,'rb')
    nth = 0
    while True:
        line = f.readline()
	
        if not line.endswith('\n'): break #Alignment file truncated or still being written
    
	if line.startswith('s'):
	    seqs.append(line.strip().split())
	
	if not line.strip() and seqs:
	    assert len(seqs) == 2
	    
	    ref_s, ref_name, ref_start, ref_size, ref_strand, ref_src_size, ref_text = seqs[0]
	    read_s, read_name, read_start, read_size, read_strand, read_src_size, read_text = seqs[1]
	    seqs = [ ]
	    
	    assert ref_strand == '+'
	    
	    forward = read_strand=='+'
	    start = int(ref_start)
	    end = start + int(ref_size) - 1
	    
	    nth += 1
	    if nth % 1000 == 0:
        	sys.stderr.write('Loading hits: %d            \r' % nth)
		sys.stderr.flush()

	    yield ref_name, read_name, read_strand=='+', start, end, ref_text.upper(), read_text.upper()
	    	
	if not line: 
	    break

def iter_hit_file_eland(filename):        
    nth = 0
    for line in open(filename, 'rb'):
        if not line.endswith('\n'): break #Alignment file truncated or still being written
	
	parts = line.rstrip().split('\t')
	
	if len(parts) < 10: continue #No unique best hit
	
	read_name, read_seq, match_type, n_U0, n_U1, n_U2, ref_filename, ref_position, direction, \
	how_Ns_interpreted = parts[:10]
	substitutions = parts[10:]
	
	forward = (direction == 'F')
	start = int(ref_position)-1
	end = start+len(read_seq)
	
	ref_seq = read_seq
	for item in substitutions:
	    base = item[-1]
	    position = int(item[:-1])-1
	    ref_seq = ref_seq[:position] + base + ref_seq[position+1:]

        if not forward:
            read_seq = sequence.reverse_complement_string(read_seq)
	    ref_seq = sequence.reverse_complement_string(ref_seq)

	nth += 1
	if nth % 1000 == 0:
            sys.stderr.write('Loading hits: %d            \r' % nth)
	    sys.stderr.flush()
	
	yield ref_filename, read_name, forward, start, end, ref_seq.upper(), read_seq.upper()


def iter_hit_file(filename):
    first_line = open(filename,'rb').readline()
    if first_line.startswith('##maf'):
        return iter_hit_file_maf(filename)
    if first_line.startswith('>'):
        return iter_hit_file_eland(filename)
    return iter_hit_file_myrialign(filename)


def read_files(argv):
    clip_start, argv = get_option_value(argv, '-s', int, 0)
    clip_end, argv = get_option_value(argv, '-e', int, 0)

    if len(argv) < 2:
        raise Bad_option('Expected at least two filenames, a reference genome and and alignment file')

    reference = sequence.sequence_file_iterator(argv[0]).next()[1] 

    #read_hits = { }
    
    hits = Hits()

    #nth = 0
    for filename in argv[1:]:
	#for line in open(filename,'rb'):
        #    if not line.endswith('\n'): continue
	#    if line.startswith('#'): continue
	
	for ref_name, name, forward, start, end, read_ali, ref_ali \
	        in iter_hit_file(filename):

	    #hit = Hit()
	    #hit.name, hit.direction, hit.n_errors, span, hit.read_ali, hit.ref_ali = line.rstrip().split()
	    i = hits.length
	    hits.resize(i+1)
	    #hits.name[i], direction, n_errors, span, hits.read_ali[i], hits.ref_ali[i] = line.rstrip().split()
	    #start, end = span.split('..')
	    #hits.start[i] = int(start)-1
	    #hits.end[i] = int(end)
	    #hits.n_errors[i] = int(n_errors)	    
	    #hits.forward[i] = (direction == 'fwd')
	    
	    hits.name[i] = name
	    hits.forward[i] = forward
	    hits.start[i] = start
	    hits.end[i] = end
	    hits.read_ali[i] = read_ali
	    hits.ref_ali[i] = ref_ali

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

	    #nth += 1
	    #if nth % 10000 == 0:
        	#sys.stderr.write('Loading hits: %d            \r' % nth)
		#sys.stderr.flush()

    hits.you_are_dirty()
    
    return reference, hits




def artplot(argv):
    try:
        only_single, argv = get_option(argv, '-u')
	prefix, argv = get_option_value(argv, '-p', lambda x:x, 'artplot')
        reference, hits = read_files(argv)
    except Bad_option, error:
        print >> sys.stderr, ''
	print >> sys.stderr, 'myr artplot [options] <reference genome> <alignments> [<alignments>...]'
	print >> sys.stderr, ''
	print >> sys.stderr, 'Alignments can be the output from "myr align", the output'
	print >> sys.stderr, 'of BLAT in "maf" format, or an ELAND results file.'
	print >> sys.stderr, ''
	print >> sys.stderr, 'Options:'
	print >> sys.stderr, ''
        print >> sys.stderr, '    -p xx - Prefix for output files, default "artplot"'
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

                in_bounds = (0 <= pos < size)
        	if a == '-':
		    if in_bounds:
		        deletions[pos] += weight
		        coverage[pos] += weight
		    pos += 1
		elif b == '-':
		    if in_bounds:
		        insertions[pos] += weight
		else:
		    if in_bounds:
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

    normalizer = numpy.maximum(1.0, coverage)

    save(prefix+'-coverage.txt', coverage)
    save(prefix+'-insertions.txt', insertions / normalizer)
    save(prefix+'-deletions.txt', deletions / normalizer)
    save(prefix+'-substitutions.txt', substitutions / normalizer)
    
    base_counts = base_counts[:,:4]
    base_total = numpy.sum(base_counts, 1)
    surprise = numpy.zeros((size,4),'float64')
    for i in xrange(4):
        good = base_counts[:,i] > 0
        surprise[good,i] = numpy.log( base_counts[good,i] / base_total[good] ) / numpy.log(0.5)
    entropy = numpy.zeros(size, 'float64')
    good = base_total > 0
    entropy[good] = numpy.sum(base_counts[good,:] * surprise[good,:], 1) / base_total[good]
    save(prefix+'-confusion.txt', entropy)


def textdump(argv):
    try:
        only_single, argv = get_option(argv, '-u')
        reference, hits = read_files(argv)
    except Bad_option, error:
        print >> sys.stderr, ''
	print >> sys.stderr, 'myr textdump [options] <reference genome> <alignments> [<alignments>...]'
	print >> sys.stderr, ''
	print >> sys.stderr, 'Alignments can be the output from "myr align", the output'
	print >> sys.stderr, 'of BLAT in "maf" format, or an ELAND results file.'
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





# =====================================================================
# (perhaps move somewhere else...)

"""

A location is <sequence id, 32 bits> <is forward?, 1 bit> <location, 31 bits>

Tables:

    Sequence:
    - name (string) - urlish, eg read:blah, cds:blahblah
    - sequence (array)

    Alignment:
    - type (string)
    
    Base_link
    - location1
    - location2
    - alignment id

"""



class Manymany:
    def __init__(self):
        self.forward = { }
        self.back = { }
    
    def create_forward(self, item):
        if item not in self.forward: 
            self.forward[item] = sets.Set()
    
    def create_back(self, item):
        if item not in self.back:
            self.back[item] = sets.Set()
    
    def create(self, item):
        self.create_forward(item)
        self.create_back(item)

    def destroy_forward(self, item):
        assert not self.forward[item]
        del self.forward[item]

    def destroy_back(self, item):
        assert not self.back[item]
        del self.back[item]
          
    def destroy(self, item):
        self.destroy_forward(item)
        self.destroy_back(item)
    
    def link(self, a, b):
        self.forward[a].add(b)
        self.back[b].add(a)
        
    def unlink(self, a, b):
        self.forward[a].remove(b)
        self.back[b].remove(a)
        

class Union:
    def __init__(self):
        self.parent = { }
    
    def create(self, item):
        if item not in self.parent:
	    self.parent[item] = item
    
    def root(self, item):
        if self.parent[item] == item:
	    return item
	self.parent[item] = self.root(self.parent[item])
	return self.parent[item]
	
    def merge_if_created(self, a,b):
        if a in self.parent and b in self.parent:
            self.parent[self.root(a)] = self.root(b)
    
    def sets(self):
        results = { }
        for item in self.parent:
	    root = self.root(item)
	    if root not in results:
	        results[root] = sets.Set()
	    results[root].add(item)
	return list(results.values())


FORWARD_MASK = numpy.uint64( (1<<31) )
POSITION_MASK = numpy.uint64( (1<<31)-1 )
ONE_UINT64 = numpy.uint64(1)
ZERO_UINT64 = numpy.uint64(0)
SEQUENCE_SHIFT = numpy.uint64(32)

def make_location(seq_id, is_forward, position):
    return (numpy.uint64(seq_id) << SEQUENCE_SHIFT) + (is_forward and FORWARD_MASK or ZERO_UINT64) + numpy.uint64(position)

def location_parts(location):
    return location>>SEQUENCE_SHIFT, bool(location&FORWARD_MASK), location&POSITION_MASK

def location_sequence(location):
    return location>>SEQUENCE_SHIFT

def location_next(location):
    """ Note: no bounds checking """
    if location & FORWARD_MASK:
        return location + ONE_UINT64
    else:
        return location - ONE_UINT64

class Sequences(Table):
    MEMBERS = (
        ('name', 'object'),
	('sequence', 'object'),
	('comment', 'object'),
    )

class Alignments(Table):
    MEMBERS = (
        ('type', 'object'),
    )

class Base_links(Table):
    MEMBERS = (
        ('location1', 'uint64'),
	('location2', 'uint64'),
	('alignment', 'uint32'),
    )



class Browser:
    def __init__(self):
        self.sequences = Sequences()
	self.name_to_sequence = { }
	self.alignments = Alignments()
	self.base_links = Base_links()
	
    def open_screen(self):
	import curses
	
	self.screen = curses.initscr()
	curses.noecho()
	curses.cbreak()
	self.screen.keypad(1)
    
    def close_screen(self):
        import curses
	
	self.screen.keypad(0)
	curses.nocbreak()
	curses.echo()
	curses.endwin()
	
    def location_move(self, location, offset):
        seq, forward, pos = location_parts(location)
        length = len(self.sequences.sequence[seq])
	if not forward:
	    offset = -offset
	pos += offset
	if not 0 <= pos < length:
	    raise Out_of_bounds()
	return make_location(seq,forward,pos)

    def location_get(self, location):
        seq, forward, pos = location_parts(location)
	result = self.sequences.sequence[seq][pos]
	if forward:
	    return result
	else:
	    return sequence.COMPLEMENT[result]

    def valid_location(self, location):
        seq, forward, pos = location_parts(location)
	return pos < len(self.sequences.sequence[seq])
	
    def add_sequence(self, name, sequence, comment=''):
        i = self.sequences.new_id()
	self.sequences.name[i] = name
	self.sequences.sequence[i] = sequence
	self.sequences.comment[i] = comment
	self.sequences.you_are_dirty()
	self.name_to_sequence[name] = i
	return i

    def add_alignment(self, type, 
                      name1, fwd1, start1, ali1, 
		      name2, fwd2, start2, ali2, reverse_counts_from_end=True):
	try:
            seq1 = self.name_to_sequence[name1]
	except KeyError:
	    raise Error('Sequence "%s" referenced by an alignment has not been loaded' % name1)
	
	try:
            seq2 = self.name_to_sequence[name2]
	except KeyError:
	    raise Error('Sequence "%s" referenced by an alignment has not been loaded' % name2)

        if reverse_counts_from_end and not fwd1:
	    start1 = len(self.sequences.sequence[seq1])-1-start1
        if reverse_counts_from_end and not fwd2:
	    start2 = len(self.sequences.sequence[seq2])-1-start2
	
	location1 = make_location(seq1, fwd1, start1)
	location2 = make_location(seq2, fwd2, start2)
	
	ali = self.alignments.new_id()
	self.alignments.type[ali] = type
	
        assert len(ali1) == len(ali2)

	i = 0
	if len(ali1):
	    while True:
		if ali1[i] != '-' and ali2[i] != '-':
	            link = self.base_links.new_id()
		    assert self.valid_location(location1), 'Bad alignment %d <%s>-%s' % (i,name1,name2)
		    assert self.valid_location(location2), 'Bad alignment %d %s-<%s>' % (i,name1,name2)
		    self.base_links.location1[link] = location1
		    self.base_links.location2[link] = location2
		    self.base_links.alignment[link] = ali

		if i == len(ali1)-1:
	            break

		if ali1[i] != '-':
	            location1 = location_next(location1)
		if ali2[i] != '-':
	            location2 = location_next(location2)

		i += 1
	
	self.alignments.you_are_dirty()
	self.base_links.you_are_dirty()
	return ali

    def load_sequences(self, filename):
        for name, seq in sequence.sequence_file_iterator(filename):
	    self.add_sequence(name, seq)

    def load_myr_hits(self, filename):
	for ref_name, name, forward, start, end, read_ali, ref_ali \
	        in iter_hit_file(filename):

	    if name not in self.name_to_sequence:
	        seq = sequence.sequence_from_string(read_ali.replace('-',''))
		if not forward:
		    seq = sequence.reverse_complement(seq)
		self.add_sequence(name, seq)		

	    seq = self.sequences.sequence[self.name_to_sequence[name]]
	    
	    self.add_alignment('myr align',
	        ref_name, True, start, ref_ali,
		name, forward, 0, read_ali)
		
	    #if forward:
	    #    read_start = 0
	    #else:
	    #    read_start = len(seq)-1
	    
	    #self.add_alignment('myr align',
	    #    ref_name, True, start, ref_ali,
	#	name, forward, read_start, read_ali)

    def load_maf(self, filename):
	seqs = [ ]
	f = open(filename,'rb')
	nth = 0
	while True:
	    line = f.readline()

	    if not line.endswith('\n'): break #Alignment file truncated or still being written

	    if line.startswith('s'):
		s, name, start, size, strand, src_size, text = line.strip().split()
		forward = strand == '+'
		start = int(start)
		seqs.append((name, forward, start, text))

	    if not line.strip() and seqs:
	        for i in xrange(len(seqs)):
		    for j in xrange(i):
		        self.add_alignment('maf',
			    seqs[i][0],seqs[i][1],seqs[i][2],seqs[i][3],
			    seqs[j][0],seqs[j][1],seqs[j][2],seqs[j][3],
			    True)
		seqs = []
		
	    if not line: 
		break

    def load_velvet_graph(self, filename):
        comments = { }
        f = open(os.path.join(filename,'stats.txt'), 'rb')
	f.readline()
	for line in f:
	    ID, lgth, n_out, n_in, long_cov, short1_cov, short1_Ocov, short2_cov, short2_Ocov = line.strip().split()
	    comments['NODE_'+ID] = 'cov=%.1f' % (float(long_cov)+float(short1_cov)+float(short2_cov))
	
        f = open(os.path.join(filename,'LastGraph'), 'rb')
        line = f.readline()
	hash_size = int( line.split()[2] )	
	tail_size = hash_size - 1
	
	while True:
	    line = f.readline()
	    if not line: break
	    parts = line.strip().split()
	    
	    if parts[0] == 'NODE':
	        node_name = 'NODE_' + parts[1]
		fwd = sequence.sequence_from_string(f.readline().strip())
		rev = sequence.sequence_from_string(f.readline().strip())
		assert len(fwd) == len(rev)
		if len(fwd) < tail_size:
		    pad = [4]*(tail_size-len(fwd))
		    fwd = numpy.concatenate((pad,fwd))
		    rev = numpy.concatenate((pad,rev))
		rev_rc = sequence.reverse_complement(rev)

		#if not numpy.alltrue(numpy.equal(fwd[:-tail_size], rev_rc[tail_size:])):
		#    print node_name
                #    print fwd[:-tail_size]
		#    print rev_rc[tail_size:]
                #    print numpy.equal(fwd[:-tail_size],rev_rc[tail_size:]).astype('int')
		#seq = numpy.concatenate((rev_rc,fwd[-tail_size:]))
		
		#self.add_sequence(node_name, seq)
		#print node_name
		
		#TODO: IUPAC codes where different
                inner_fwd = fwd[:-tail_size]
		inner_rev = rev_rc[tail_size:]		
		self.add_sequence(
		  node_name, 
		  numpy.concatenate((
		    rev_rc[:tail_size],
		    numpy.where(numpy.equal(inner_fwd, inner_rev),
		                inner_fwd,
				4),
		    fwd[-tail_size:]
		  )),
		  comments[node_name])
		
		#self.add_sequence(node_name+'_fwd', fwd)
		#self.add_sequence(node_name+'_rev', rev)
		#self.add_alignment('velvet_contig_pair',
		#    node_name+'_fwd', True,  0,
		#    sequence.string_from_sequence(fwd[:-tail_size]),
		#    node_name+'_rev', False, len(rev_rc)-tail_size-1,
		#    sequence.string_from_sequence(rev_rc[tail_size:]) )
	
	    if parts[0] == 'ARC':
	        node_from = int(parts[1])
		name_from = 'NODE_%d' % abs(node_from)
		fwd_from = node_from >= 0

	        node_to = int(parts[2])
		name_to = 'NODE_%d' % abs(node_to)
		fwd_to = node_to >= 0
		
		len_from = len(self.sequences.sequence[ self.name_to_sequence[name_from] ])
		
		self.add_alignment('velvet_arc',
		     name_from, fwd_from, len_from-tail_size, 'X'*tail_size,
		     name_to, fwd_to, 0, 'X'*tail_size)
	    
	        #node_from = int(parts[1])		
		#if node_from < 0:
		#    node_from = 'NODE_%s_rev' % -node_from
		#else:
		#    node_from = 'NODE_%s_fwd' % node_from
                #node_from_id = self.name_to_sequence[node_from]		
		
		#node_to = -int(parts[2])
		#if node_to < 0:
		#    node_to = 'NODE_%s_rev' % -node_to
		#else:
		#    node_to = 'NODE_%s_fwd' % node_to
                #node_to_id = self.name_to_sequence[node_to]
		
		#len_from = len(self.sequences.sequence[node_from_id])
		#len_to = len(self.sequences.sequence[node_to_id])
		#self.add_alignment('velvet_arc',
		#    node_from, True,  len_from-tail_size, 'X'*tail_size,
		#    node_to,   False, len_to-1, 'X'*tail_size)
		
		
		

    def show(self, cursor, distance_cutoff):
	positions = { }
	todo = [ ]
	heapq.heappush(todo, (0, 0, cursor))
	def add_todo(location, distance, position):
	    if distance > distance_cutoff:
	        raise Out_of_bounds()
	    if location in positions: 
	        return
	    assert self.valid_location(location)
	    heapq.heappush(todo, (distance, position, location))
	
	#dag = Dag()
	dag = { }
	def dag_link(a,b):
	    if a not in dag: dag[a] = [ ]
	    if b not in dag: dag[b] = [ ]
	    dag[a].append(b)
	contigua = Union()
	
	while todo:
	    distance, position, location = heapq.heappop(todo)
	    if location in positions: continue
	    positions[location] = position
	    
	    #dag.get_keyset(location)
	    if location not in dag: dag[location] = [ ]
	    
	    contigua.create(location)
	    
	    #flipped_location = location ^ FORWARD_MASK
	    #add_todo(flipped_location, distance)
	    #dag.merge_keys(flipped_location, location)

	    try:
	        linked_location = self.location_move(location,1)
	        contigua.merge_if_created(location, linked_location)
	        add_todo(linked_location, distance+1, position+1)
		dag_link(location, linked_location)
	    except Out_of_bounds: pass

	    try:
	        linked_location = self.location_move(location,-1)
	        contigua.merge_if_created(location, linked_location)
	        add_todo(linked_location, distance+1, position-1)
		dag_link(linked_location, location)
	    except Out_of_bounds: pass
	    
	    def merge(linked_location):
	        try:
		    add_todo(linked_location, distance, position)
		    dag_link(location, linked_location)
		    dag_link(linked_location, location)
		except Out_of_bounds: pass
	    
	    for i in self.base_links.find_all('location1',location):
	        merge(self.base_links.location2[i])
	    for i in self.base_links.find_all('location2',location):
	        merge(self.base_links.location1[i])

	    for i in self.base_links.find_all('location1',location^FORWARD_MASK):
	        merge(self.base_links.location2[i]^FORWARD_MASK)
	    for i in self.base_links.find_all('location2',location^FORWARD_MASK):
	        merge(self.base_links.location1[i]^FORWARD_MASK)


	
	class Contig: pass
	contigs = [ ]
	for item in contigua.sets():
	    sample = iter(item).next()
	    seq = location_sequence( sample )
	    forward = (sample & FORWARD_MASK) != 0
	    contig = Contig()
	    contigs.append(contig)
	    contig.seq = seq
	    contig.name = self.sequences.name[seq]
	    contig.forward = forward
	    contig.sort_key = (contig.name, not forward)
	    contig.locations = item
	    
	contigs.sort(lambda a,b: cmp(a.sort_key, b.sort_key))
		
	#print contigua
	
	#for contig in contigua:
	#    item = iter(contig).next()
	#    seq = location_sequence( item )
	#    forward = (item & FORWARD_MASK) != 0
	#    print self.sequences.name[seq], forward
	
	#order = dag.sort(positions)
	def priority(component):
	    return float(sum([ positions[item] for item in component ])) / len(component)
	order = sort.compact_robust_topological_sort(dag, priority)
	
	table = [ ]
	column_width = [ ]
	
	for x, locations in enumerate(order):
	    column = [ ]
	    table.append(column)
	    for y, contig in enumerate(contigs):
	        relevant = [ location for location in locations
		             if location in contig.locations ]
		relevant.sort()
		if relevant and not (relevant[0]&FORWARD_MASK):
		    relevant = relevant[::-1]
		if self.cursor in relevant:
		    cursor_y = y
		    cursor_x = numpy.sum(column_width) + relevant.index(self.cursor)
		    cursor_column = column
	        column.append( relevant )
		#sequence.string_from_sequence( [ self.location_get(location)
		#                          for location in relevant ] ) )
	    
	    column_width.append(max([ len(item) for item in column ]))


        self.screen.clear()
	
	maxy, maxx = self.screen.getmaxyx()
	offset_y = int( maxy//2-cursor_y )
	offset_x = int( maxx//2-cursor_x )
	def addstr(y,x,string):
	    if y < 0 or y >= maxy: return
	    while string and x < 0:
	        string = string[1:]
		x += 1
	    if x+len(string) > maxx:
	        string = string[:max(0,maxx-x)]
	    if not string: return
	    #try:
	    self.screen.addstr(y,x,string)
	    #except:
	    #    raise repr((y,x,string,maxy,maxx))
	
        for y in xrange(len(contigs)):
	
	    #item = iter(contigua[y]).next()
	    #seq = location_sequence( item )
	    #forward = (item & FORWARD_MASK) != 0
	    #sys.stdout.write('% 20s %d  ' % (self.sequences.name[seq],forward))
	    
	    scr_x = 0
            for x in xrange(len(order)):
	        item = table[x][y]
		#item += ' '*(column_width[x]-len(item))
	        #sys.stdout.write(item)
		
		string = sequence.string_from_sequence( [ self.location_get(location)
		                                          for location in item ] )
		
		addstr(y+offset_y,scr_x+offset_x,string)
		
		scr_x += column_width[x]
		
	    #sys.stdout.write('\n')
	    
	    info = contigs[y].name
	    if self.sequences.comment[contigs[y].seq]:
	        info += ' ' + self.sequences.comment[contigs[y].seq]
	    if contigs[y].forward:
	        info += ' >>> '
	    else:
	        info += ' <<< '
	    addstr(y+offset_y,max(0,-len(info)-1+offset_x),info)
	    
	cursor_seq, cursor_fwd, cursor_pos = location_parts(cursor)
	addstr(1,1, '%s @ %d' % (self.sequences.name[cursor_seq], cursor_pos))
	    
	self.screen.move(cursor_y+offset_y,cursor_x+offset_x)
	self.screen.refresh()
	
	return cursor_column, cursor_y
	
    def browse(self, initial=None):
        import curses
    
        if self.sequences.length == 0:
	    raise Error('No sequences to browse')

        if initial is not None:
	    initial_id = self.name_to_sequence[initial]
	else:
	    initial_id = 0
        self.cursor = make_location(initial_id,True,0)
	
	radius = 60
        cursor_column, cursor_y = self.show(self.cursor, radius)
	while True:
	    self.screen.nodelay(1)
	    key = self.screen.getch()
	    self.screen.nodelay(0)
	
	    if key == -1:
	        cursor_column, cursor_y = self.show(self.cursor, radius)
	        key = self.screen.getch()
	    
	    if key == curses.KEY_LEFT:
	        try:
		    self.cursor = self.location_move(self.cursor, -1)
		except Out_of_bounds:
		    pass
	    elif key == curses.KEY_RIGHT:
	        try:
		    self.cursor = self.location_move(self.cursor, 1)
		except Out_of_bounds:
		    pass
	    elif key == curses.KEY_UP:
	        while True:
		    cursor_y -= 1
		    if cursor_y < 0 or cursor_column[cursor_y]: break
		if cursor_y >= 0: self.cursor = cursor_column[cursor_y][0]		    
	    elif key == curses.KEY_DOWN:
	        while True:
		    cursor_y += 1
		    if cursor_y >= len(cursor_column) or cursor_column[cursor_y]: break
		if cursor_y < len(cursor_column): self.cursor = cursor_column[cursor_y][0]		    
	    elif key == 10:
	        self.cursor = self.cursor ^ FORWARD_MASK
	    elif key == curses.KEY_HOME:
	        seq_id, fwd, pos = location_parts(self.cursor)
		if fwd:
		    pos = 0
		else:
		    pos = len(self.sequences.sequence[seq_id])-1
		self.cursor = make_location(seq_id, fwd, pos)
	    elif key == curses.KEY_END:
	        seq_id, fwd, pos = location_parts(self.cursor)
		if not fwd:
		    pos = 0
		else:
		    pos = len(self.sequences.sequence[seq_id])-1
		self.cursor = make_location(seq_id, fwd, pos)
	    elif key == 27:
	        break
	    
	
        #for i in xrange(0,300,1):	    
        #    self.cursor = make_location(0,True,i)
	#    self.show(self.cursor, 10)
	    #import time; time.sleep(0.1)

BROWSE_USAGE = """\

myr browse [sequence files] -aligns [alignment files]

myr browse -velvet [velvet dir]

"""

def browse(argv):
    if not argv:
        sys.stderr.write(BROWSE_USAGE)
	return 1

    browser = Browser()
    
    modes = ['-seqs', '-aligns', '-velvet','-maf','-initial']
    mode = '-seqs'
    initial = None
    for item in argv:
        if item in modes:
	    mode = item
	elif mode == '-seqs':
	    browser.load_sequences(item)
	elif mode == '-aligns':
	    browser.load_myr_hits(item)
	elif mode == '-velvet':
	    browser.load_velvet_graph(item)
	elif mode == '-maf':
	    browser.load_maf(item)
	elif mode == '-initial':
	    assert initial is None, 'More than one initial sequence name given'
	    initial = item

    browser.open_screen()
    try:
	browser.browse(initial)    
    finally:
        browser.close_screen()

    return 0


