
import numpy, sys
from numpy import random

import sequence

def main(argv):
    if len(argv) != 2:
        print 
        print 'myr shred'
        print
        print 'Generate fake Illumina reads.'
        print 'Not guaranteed to be sanely calibrated, for testing only.'
        print
        print 'Usage:'
        print
        print '    myr shred <number of reads> <sequence.fna>'
        print
        return 1
         
    how_many = int(argv[0])
    seq = sequence.sequence_file_iterator(argv[1]).next()[1]

    READ_SIZE = 33
    error_p = numpy.array(
          [ 0.00912327,  0.00930828,  0.00929492,  0.00928049,  0.0093261 ,
            0.00928905,  0.00938066,  0.00936397,  0.00939301,  0.00947136,
            0.00952966,  0.00956763,  0.01073044,  0.01091972,  0.01121085,
            0.01159389,  0.01200634,  0.01233303,  0.01271543,  0.01334389,
            0.01349712,  0.01412138,  0.01462227,  0.01720922,  0.01617627,
            0.01671721,  0.01795653,  0.01904574,  0.02032015,  0.0220367 ,
            0.02354595,  0.02560759,  0.03480737])
    
    for i in xrange(how_many):
        print '>read%d' % i
    
        pos = random.randint(len(seq)-READ_SIZE+1)
        read = seq[pos:pos+READ_SIZE]
        if random.randint(2): read = sequence.reverse_complement(read)
        
        read = read.copy()
        mutations = random.random(READ_SIZE) < error_p
        read[mutations] = (
            read[mutations] + 
    	random.randint(1,4,size=numpy.sum(mutations)).astype('uint8')
        ) % 4
    
        print sequence.string_from_sequence(read)
    
