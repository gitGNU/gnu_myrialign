
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

import sys, numpy, os.path, sets, heapq

import sequence

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
	('n_errors', 'int32'),
    )   


def iter_hit_file(filename):
    ref_name = None
    nth = 0
    for line in open(filename, 'rb'):
        if line.startswith('#'):
            if line.startswith('#Reference:'):
	        ref_name = line.split()[1]
	    continue
	
	name, direction, n_errors, span, read_ali, ref_ali = line.rstrip().split()
	start, end = span.split('..')
	start = int(start)-1
	end = int(end)
	n_errors = int(n_errors)	    
	forward = (direction == 'fwd')
	
	nth += 1
	if nth % 1000 == 0:
            sys.stderr.write('Loading hits: %d            \r' % nth)
	    sys.stderr.flush()
	
	yield ref_name, name, forward, n_errors, start, end, read_ali, ref_ali
	

def read_files(argv):
    clip_start, argv = get_option_value(argv, '-s', int, 0)
    clip_end, argv = get_option_value(argv, '-e', int, 0)

    if len(argv) < 2:
        raise Bad_option('Expected at least two filenames, a reference genome and output from myr align')

    reference = sequence.sequence_file_iterator(argv[0]).next()[1] 

    #read_hits = { }
    
    hits = Hits()

    #nth = 0
    for filename in argv[1:]:
	#for line in open(filename,'rb'):
        #    if not line.endswith('\n'): continue
	#    if line.startswith('#'): continue
	
	for ref_name, name, forward, n_errors, start, end, read_ali, ref_ali \
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
	    hits.n_errors[i] = n_errors
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
        

class Dag(Manymany):
    def __init__(self):
        Manymany.__init__(self)
        self.key_keyset = Manymany()
    
    def has_key(self, key):
        return key in self.key_keyset.forward
    
    def get_keyset(self, key):        
        if key not in self.key_keyset.forward:
            keyset = (key,)
            self.key_keyset.create_forward(key)
            self.key_keyset.create_back(keyset)
            self.key_keyset.link(key, keyset)
            self.create(keyset)
        
        return iter(self.key_keyset.forward[key]).next() #Should use a manyone, really
	
    def get_all_keysets(self):
        return self.key_keyset.back.keys()
    
    def _find_betweeners(self, keyset, visited, betweeners):
        if keyset in betweeners:
            return True
            
        if keyset in visited:
            return False
                    
        visited.add(keyset)
        
        any = False
        
        for next in self.forward[keyset]:
            any = any or self._find_betweeners(next,visited,betweeners)
        
        if any:
            betweeners.add(keyset)
            
        return any
    
    def link_keys(self, a, b):
        a = self.get_keyset(a)
        b = self.get_keyset(b)
        
        visited = sets.Set()
        betweeners = sets.Set((a,))
        if self._find_betweeners(b, visited, betweeners):
            self.merge_keysets(betweeners)
        else:
            self.link(a,b)
    
    def merge_keys(self, a, b):
        #a = self.get_keyset(a)
        #b = self.get_keyset(b)
	#if a != b:
        #    self.merge_keysets((a,b))
	
	self.link_keys(a,b)
	self.link_keys(b,a)
        
    def merge_keysets(self, keysets):
        if len(keysets) <= 1:
	    return
	    
        new_keyset = sum(keysets, ())
        self.create(new_keyset)
        self.key_keyset.create_back(new_keyset)
        
        for keyset in keysets:
            for old_keyset in self.forward[keyset].copy():
                self.unlink(keyset, old_keyset)
                if old_keyset not in keysets:
                    self.link(new_keyset, old_keyset)
            for old_keyset in self.back[keyset].copy():
                self.unlink(old_keyset, keyset)
                if old_keyset not in keysets:
                    self.link(old_keyset, new_keyset)

            for key in keyset:
                self.key_keyset.unlink(key, keyset)
                self.key_keyset.link(key, new_keyset)
            
        for keyset in keysets:
            self.key_keyset.destroy_back(keyset)  
            self.destroy(keyset)

    def sort(self):
        ready = [ ]
        counters = { }
        for keyset in self.get_all_keysets():
            counters[keyset] = len(self.back[keyset])
            if not counters[keyset]:
                ready.append(keyset)
        
        result = [ ]   
        while ready:
            item = ready.pop(0) #TODO: clever ordering: furthest before cursor, nearest after cursor
            result.append(item)
            for keyset in self.forward[item]:
                counters[keyset] -= 1
		assert counters[keyset] >= 0
                if not counters[keyset]:
                    ready.append(keyset)
        return result

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
	import curses
	
        self.sequences = Sequences()
	self.name_to_sequence = { }
	self.alignments = Alignments()
	self.base_links = Base_links()
	
	self.screen = curses.initscr()
	curses.noecho()
	curses.cbreak()
	self.screen.keypad(1)
    
    def close(self):
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
	
    def add_sequence(self, name, sequence):
        i = self.sequences.new_id()
	self.sequences.name[i] = name
	self.sequences.sequence[i] = sequence
	self.sequences.you_are_dirty()
	self.name_to_sequence[name] = i
	return i

    def add_alignment(self, type, 
                      name1, fwd1, start1, ali1, 
		      name2, fwd2, start2, ali2):
	try:
            seq1 = self.name_to_sequence[name1]
	except KeyError:
	    raise Error('Sequence "%s" referenced by an alignment has not been loaded' % name1)
	
	try:
            seq2 = self.name_to_sequence[name2]
	except KeyError:
	    raise Error('Sequence "%s" referenced by an alignment has not been loaded' % name2)
	
	location1 = make_location(seq1, fwd1, start1)
	location2 = make_location(seq2, fwd2, start2)
	
	ali = self.alignments.new_id()
	self.alignments.type[ali] = type
	
	for i in xrange(len(ali1)):
	    if ali1[i] != '-' and ali2[i] != '-':
	        link = self.base_links.new_id()
		self.base_links.location1[link] = location1
		self.base_links.location2[link] = location2
		self.base_links.alignment[link] = ali
	    
	    if ali1[i] != '-':
	        location1 = location_next(location1)
	    if ali2[i] != '-':
	        location2 = location_next(location2)
	
	self.alignments.you_are_dirty()
	self.base_links.you_are_dirty()
	return ali

    def load_sequences(self, filename):
        for name, seq in sequence.sequence_file_iterator(filename):
	    self.add_sequence(name, seq)

    def load_myr_hits(self, filename):
	for ref_name, name, forward, n_errors, start, end, read_ali, ref_ali \
	        in iter_hit_file(filename):

	    if name not in self.name_to_sequence:
	        seq = sequence.sequence_from_string(read_ali.replace('-',''))
		if not forward:
		    seq = sequence.reverse_complement(seq)
		self.add_sequence(name, seq)		

	    seq = self.sequences.sequence[self.name_to_sequence[name]]
	    
	    if forward:
	        read_start = 0
	    else:
	        read_start = len(seq)-1
	    
	    self.add_alignment('myr align',
	        ref_name, True, start, ref_ali,
		name, forward, read_start, read_ali)

    def show(self, cursor, distance_cutoff):
	done = sets.Set()
	todo = [ ]
	heapq.heappush(todo, (0, cursor))
	def add_todo(location, distance):
	    if distance > distance_cutoff:
	        raise Out_of_bounds()
	    if location in done: 
	        return
	    heapq.heappush(todo, (distance, location))
	
	dag = Dag()
	contigua = Union()
	
	while todo:
	    distance, location = heapq.heappop(todo)
	    if location in done: continue
	    done.add(location)
	    
	    dag.get_keyset(location)
	    
	    contigua.create(location)
	    
	    #flipped_location = location ^ FORWARD_MASK
	    #add_todo(flipped_location, distance)
	    #dag.merge_keys(flipped_location, location)

	    try:
	        linked_location = self.location_move(location,1)
	        contigua.merge_if_created(location, linked_location)
	        add_todo(linked_location, distance+1)
		dag.link_keys(location, linked_location)
	    except Out_of_bounds: pass

	    try:
	        linked_location = self.location_move(location,-1)
	        contigua.merge_if_created(location, linked_location)
	        add_todo(linked_location, distance+1)
		dag.link_keys(linked_location, location)
	    except Out_of_bounds: pass
	    
	    def merge(linked_location):
	        try:
		    add_todo(linked_location, distance+1)
		    dag.merge_keys(location, linked_location)
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
	
	order = dag.sort()
	
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
	    
	    column_width.append(max([ len(item) for item in column ]) + 1)


        self.screen.clear()
	
	maxy, maxx = self.screen.getmaxyx()
	offset_y = maxy//2-cursor_y
	offset_x = maxx//2-cursor_x
	
        for y in xrange(len(contigs)):
	    info = contigs[y].name
	    if contigs[y].forward:
	        info += ' >>>'
	    else:
	        info += ' <<<'
	    self.screen.addstr(y+offset_y,-len(info)-1+offset_x,info)
	
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
		
		self.screen.addstr(y+offset_y,scr_x+offset_x,string)
		
		scr_x += column_width[x]
		
	    #sys.stdout.write('\n')
	    
	self.screen.move(cursor_y+offset_y,cursor_x+offset_x)
	self.screen.refresh()
	
	return cursor_column, cursor_y
	
    def browse(self):
        import curses
    
        if self.sequences.length == 0:
	    raise Error('No sequences to browse')

        self.cursor = make_location(0,True,0)
	
        cursor_column, cursor_y = self.show(self.cursor, 20)
	while True:
	    self.screen.nodelay(1)
	    key = self.screen.getch()
	    self.screen.nodelay(0)
	
	    if key == -1:
	        cursor_column, cursor_y = self.show(self.cursor, 20)
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
	    elif key == 27:
	        break
	    
	
        #for i in xrange(0,300,1):	    
        #    self.cursor = make_location(0,True,i)
	#    self.show(self.cursor, 10)
	    #import time; time.sleep(0.1)

BROWSE_USAGE = """\

myr browse [sequence files] -aligns [alignment files]

"""

def browse(argv):
    if not argv:
        sys.stderr.write(BROWSE_USAGE)
	return 1

    browser = Browser()
    try:
    
        modes = ['-seqs', '-aligns']
        mode = '-seqs'
        for item in argv:
            if item in modes:
	        mode = item
	    elif mode == '-seqs':
	        browser.load_sequences(item)
	    elif mode == '-aligns':
		browser.load_myr_hits(item)

	browser.browse()
    
    finally:
        browser.close()

    return 0


