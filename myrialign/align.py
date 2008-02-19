
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
    
    Exploits bit-parallelism, as well as conventional multi-processor 
    parallelism to align short reads to a reference.
    
    Uses Cell processor SPUs if available.

"""

import numpy, random, time, sys, os, string, select, struct, fcntl

import spu, children, sequence

def how_many_cpus():
    """Detects the number of effective CPUs in the system,
    
       Function nicked from Parallel Python."""
    #for Linux, Unix and MacOS
    if hasattr(os, "sysconf"):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"): 
            #Linux and Unix
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: 
            #MacOS X
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    #for Windows
    if os.environ.has_key("NUMBER_OF_PROCESSORS"):
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
        if ncpus > 0:
            return ncpus
    #return the default value
    return 1

def is_cell():
    if os.path.exists('/proc/cpuinfo'):
        for line in open('/proc/cpuinfo','rU'):
            if line.startswith('cpu') and \
               'Cell Broadband Engine' in line:
                return True
    return False

# TODO: Be cleverer about Cell SPU count
CELL_PROCESSOR = is_cell()
if CELL_PROCESSOR:
    PROCESSES = 7 #PS3: 6 SPUs and one thread to get ready
else:
    PROCESSES = how_many_cpus()
    


# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
#                             Bit vector utilities
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================

TYPE = 'uint64'
BITS = 64
#TYPE = 'uint32'
#BITS = 32
TRUE = numpy.array(-1, TYPE)

#Little endian
#BIT = numpy.array(1,TYPE) << numpy.arange(BITS, dtype=TYPE)

#Big endian
BIT = numpy.array(1,TYPE) << numpy.arange(BITS-1,-1,-1, dtype=TYPE)

def collapse(array, block_size=BITS):
    shape = array.shape
    ndim = len(shape)
    outshape = shape[:ndim-1] + ((shape[-1]+block_size-1)//block_size*block_size//BITS ,)
    output = numpy.zeros(outshape, TYPE)
    for i in xrange(BITS):
        this_bit = array[...,i::BITS]
        output[...,:this_bit.shape[-1]][this_bit] |= BIT[i]
    return output

def expand(array):
    shape = array.shape
    ndim = len(shape)
    outshape = shape[:ndim-1] + (shape[-1]*BITS,)
    output = numpy.zeros(outshape, 'bool')
    for i in xrange(BITS):
        output[...,i::BITS] = (array & BIT[i]) != 0
    return output



# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
#                             Sequence utilities
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================

#_complement = string.maketrans('ATGC','TACG')
#def reverse_complement(seq):
#    return string.translate(seq,_complement)[::-1]
#
#_nuc_trans = string.maketrans('ACGT',chr(0)+chr(1)+chr(2)+chr(3))
#def sequence_nucs(sequence):
#    return numpy.fromstring(
#        string.translate(sequence, _nuc_trans), 'uint8')
#    #sequence = numpy.fromstring(sequence.upper(), 'uint8')
#    #result = numpy.zeros(len(sequence), 'uint8')
#    #result[sequence == ord('A')] = 0
#    #result[sequence == ord('C')] = 1
#    #result[sequence == ord('G')] = 2
#    #result[sequence == ord('T')] = 3
#    #return result
#
#def sequence_nuc_str(sequence):
#    return string.translate(sequence, _nuc_trans)

#def sequence_nucmatch(sequence): # [ nuc, position ]
#    sequence = numpy.fromstring(sequence.upper(), 'uint8')
#    return numpy.array([
#        sequence == ord('A'),
#        sequence == ord('C'),
#        sequence == ord('G'),
#        sequence == ord('T'),
#    ])

def sequence_nucmatch(sequence): # [ nuc, position ]
    #TODO: handle Ns with greater memory efficiency
    return numpy.array([
        sequence == 0,
        sequence == 1,
        sequence == 2,
        sequence == 3,
        numpy.zeros(len(sequence), 'bool')
    ])

def align(seq1, seq2):
    """ Produce an alignment (for once we have found a hit).  
        Start point is zero in both seqs.
        End point may be anywhere in seq2, must be end of seq1. """
    len1 = len(seq1)
    len2 = len(seq2)
    scores = numpy.zeros((len1+1,len2+1),'int')
    scores[0,:] = numpy.arange(len2+1)
    scores[:,0] = numpy.arange(len1+1)
    
    #match_score = 1 - sequence.EQUAL(numpy.fromstring(seq1,'uint8')[:,None],
    #                                 numpy.fromstring(seq2,'uint8')[None,:])
    
    for i in xrange(1,len1+1):        
        for j in xrange(1,len2+1):
            #if seq1[i-1] == seq2[j-1]: match_score = 0
            #else: match_score = 1
            
            scores[i,j] = min(
                scores[i-1,j-1] + sequence.NOTEQUAL[seq1[i-1],seq2[j-1]],
                scores[i-1,j] + 1,
                scores[i,j-1] + 1
            )
    
    end2 = numpy.argmin(scores[len1,:])
    
    str_seq1 = sequence.string_from_sequence(seq1)
    str_seq2 = sequence.string_from_sequence(seq2)
    
    pos1 = len1
    pos2 = end2
    ali1 = [ ]
    ali2 = [ ]
    while True:
        if pos1 and pos2:
            step = scores[pos1-1,pos2-1]
            del1 = scores[pos1-1,pos2]
            del2 = scores[pos1,pos2-1]
            if step <= del1 and step <= del2:
                ali1.append(str_seq1[pos1-1])
                ali2.append(str_seq2[pos2-1])
                pos1 -= 1
                pos2 -= 1
            elif del1 <= del2:
                ali1.append(str_seq1[pos1-1])
                ali2.append('-')
                pos1 -= 1
            else:
                ali1.append('-')
                ali2.append(str_seq2[pos2-1])
                pos2 -= 1
        elif pos1:
            ali1.append(str_seq1[:pos1])
            ali2.append('-'*pos1)
            break
        else:
            ali1.append('-'*pos2)
            ali2.append(str_seq2[:pos2])
            break
    
    return ''.join(ali1[::-1]), ''.join(ali2[::-1]), end2, scores[len1,end2]


# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
#                       Bit-parallel alignment
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================

def observe(matchin,matchout, nucmatches):
    n_errors, n_positions  = matchin.shape[:2]
    
    matchout[0,:] = nucmatches[:]
    matchout[0,1:] &= matchin[0,:-1]
        
    for i in xrange(1,n_errors):
        matchout[i,i:] = nucmatches[i:]
        matchout[i,i+1:] &= matchin[i,i:-1]

        #Mismatch
        matchout[i,i:] |= matchin[i-1,i-1:-1]        

        # Deletion in read
        matchout[i,i:] |= matchin[i-1,i:]
        
        # Deletion in reference
        matchout[i,i:] |= matchout[i-1,i-1:-1]


def handle_hit(reference, ref_pos, read, read_name, n_errors, callback):
    ref_start = max(0, ref_pos - (len(read)-1) - n_errors)
    ref_scrap = reference[ref_start:ref_pos+1]
    ali_read, ali_scrap, scrap_start, ali_errors = \
        align(read[::-1], ref_scrap[::-1])
    ali_read = ali_read[::-1]
    ali_scrap = ali_scrap[::-1]
    ref_start = ref_pos+1 - scrap_start
    
    assert n_errors == ali_errors

    callback('%s %d %d..%d %s %s' % (read_name, n_errors, ref_start+1, ref_pos, ali_read, ali_scrap))
    

def search_cpu(reference, reads, read_names, maxerror, callback):
    # Reads *must* all be the same length
    readlen = len(reads[0])

    nucmatch = numpy.transpose(
        [ sequence_nucmatch(read) for read in reads ],
        (1,2,0) )
    nucmatch = collapse(nucmatch)

    match_in = numpy.zeros((maxerror+1, readlen, len(reads)), 'bool')
    for i in xrange(maxerror):
        match_in[i+1:,i,:] = True
    match_in = collapse(match_in)
    match_out = match_in.copy()
    
    for ref_pos, nuc in enumerate(reference):
        observe(match_in,match_out, nucmatch[nuc])
    
        hits = match_out[maxerror,readlen-1]
        if numpy.any(hits):
            for i in xrange(len(hits)):
                if hits[i]:
                    for j in xrange(BITS):
                        if not hits[i] & BIT[j]: continue
                        read_no = i*BITS+j
                        
                        for n_errors in xrange(maxerror+1):
                            if match_out[n_errors,readlen-1,i] & BIT[j]:
                                break

                        #Superior prior match?
                        if n_errors and match_in[n_errors-1,readlen-1,i] & BIT[j]:
                            continue
                        
                        #Equivalent future match?
                        if n_errors and match_out[n_errors-1,readlen-2,i] & BIT[j]:
                            continue

                        handle_hit(reference, ref_pos, reads[read_no], read_names[read_no], n_errors, callback)

        match_out, match_in = match_in, match_out


def search_spu(reference, reads, read_names, maxerror, callback):
    # Reads *must* all be the same length
    readlen = len(reads[0])
    
    nucmatch = numpy.transpose(
        [ sequence_nucmatch(read) for read in reads ],
        (1,2,0) )
    nucmatch = collapse(nucmatch, 128)
    n_vecs = nucmatch.shape[2] * BITS // 128

    spu_filename = spu.get_matcher(maxerror+1,readlen,n_vecs)

    #child_stdin, child_stdout = os.popen2('elfspe %s' % spu_filename); #Hmmm
    child = children.Child(['elfspe', spu_filename])
                
    child.write(nucmatch.tostring())
    child.write(reference.tostring())
    child.close_stdin()
    
    while True:
        children.wait([child])
        
        hit = child.read(12)
        if not hit: break
        
        hit_ref_pos, hit_read_no, hit_n_error = struct.unpack('lll', hit)
        handle_hit(reference, hit_ref_pos, reads[hit_read_no], read_names[hit_read_no], hit_n_error, callback)
    



# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================
#                             Main
# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================

def child(argv):
    if CELL_PROCESSOR:
        search_func = search_spu
    else:
        search_func = search_cpu

    while True:
        try:
            message, value = children.receive()
        except EOFError:
            break
        
        if message == 'align':
            reads, read_names, maxerror = value
            search_func(reference, reads, read_names, maxerror, 
                        lambda hit: children.send(('hit',hit)) )
            children.send(('done', len(reads)))
        elif message == 'ref':
            reference = value

    return 0


def main(argv):
    if len(argv) < 3:
        print
        print 'myr align'
        print
        print 'Align short reads to a reference genome.'
        print
        print 'Usage:'
        print
        print '  myr align <maximum errors> <reference genome file> <read file> [<read file>...]'
        print
        print 'Files can be in FASTA or ELAND format.'
        print
        print 'Each insertion, deletions, or subsitution counts as one error. The whole read (not'
	print 'just part of it) must align to the reference with less than the specified maximum'
	print 'errors in order to produce a hit.'
        print
        return 1

    if CELL_PROCESSOR:
        CHUNK = 2048
    else:
        CHUNK = 8192
    print >> sys.stderr, 'Myralign'
    
    if CELL_PROCESSOR:
        print >> sys.stderr, 'Cell processor detected'
    else:
        print >> sys.stderr, 'Cell processor not detected'
    
    print >> sys.stderr, 'Aligning batches of', CHUNK//2, 'reads per process,', PROCESSES, 'processes.'
    
    maxerror = int(argv[0])
    
    waiting = [ children.Self_child() for i in xrange(PROCESSES) ]
    running = [ ]
    
    t1 = time.time()
    total_alignments = [0]
    
    def handle_events():
        for child in children.wait(running):
            message, value = child.receive()
            if message == 'done':
                running.remove(child)
                waiting.append(child)
                
                dt = time.time() - t1
                total_alignments[0] += value//2 # Forwards + backwards == 1 alignment
                print >> sys.stderr, '%d alignments in %.2f seconds, %.4f per alignment' % (total_alignments[0], dt, dt/total_alignments[0])
            else:
                print value
    
    for ref_name, ref_seq in sequence.sequence_file_iterator(argv[1]):
        print '#', ref_name
        
        for child in waiting:
            child.send(('ref', ref_seq))
        
        # Collect reads of the same length,
        # and do them in batches
        buckets = { } # length -> [ [name], [seq] ]
        def do_bucket(length):
            read_names, read_seqs = buckets[length]
            del buckets[length]
        
            while not waiting: 
                handle_events()
        
            child = waiting.pop()
            child.send(('align', (read_seqs, read_names, maxerror)))
            running.append(child)
        
        for read_name, read_seq in sequence.sequence_files_iterator(argv[2:]):
            length = len(read_seq)
            if length not in buckets:
                buckets[length] = ( [], [] )
            buckets[length][0].append(read_name + ' fwd')
            buckets[length][1].append(read_seq)
            buckets[length][0].append(read_name + ' rev')
            buckets[length][1].append(sequence.reverse_complement(read_seq))
            
            if len(buckets[length][0]) >= CHUNK:
                do_bucket(length)
        
        for length in list(buckets):
            do_bucket(length)
        
        while running: 
            handle_events()
    
    for child in waiting:
        child.close()
    
    return 0

