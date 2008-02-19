
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

import numpy

class Parse_error(Exception): pass


SEQ_STR = numpy.empty(256,'uint8')
SEQ_STR[:] = 4 # N
SEQ_STR[ord('A')] = 0
SEQ_STR[ord('T')] = 1
SEQ_STR[ord('C')] = 2
SEQ_STR[ord('G')] = 3
SEQ_STR[ord('a')] = 0
SEQ_STR[ord('t')] = 1
SEQ_STR[ord('c')] = 2
SEQ_STR[ord('g')] = 3

COMPLEMENT = numpy.array([ 1,0,3,2,4 ], 'uint8')

STR_SEQ = numpy.array([ ord('A'),ord('T'),ord('C'),ord('G'),ord('N') ], 'uint8')

EQUAL = numpy.array(
 [ [1,0,0,0,0],
   [0,1,0,0,0],
   [0,0,1,0,0],
   [0,0,0,1,0],
   [0,0,0,0,0] ],'bool')

NOTEQUAL = ~EQUAL

def sequence_from_string(string):
    return SEQ_STR[ numpy.fromstring(string, 'uint8') ]

def string_from_sequence(seq):
    return STR_SEQ[ seq ].tostring()

def reverse_complement(seq):
    return COMPLEMENT[ seq[::-1] ] 


def fasta_iterator(filename):
    cur_seq_name = None    
    cur_seq = [ ]
    f = open(filename,'rU')
    while True:
        line = f.readline()

        # Sequence data?        
        if line and not line.startswith('>'):
            assert cur_seq_name is not None
            cur_seq.append(line.strip())
            continue
        
        # Ok, not sequence data
	# Must be the start of a new sequence or end of file
		
	if cur_seq_name is None and cur_seq:
	    raise Parse_error()
	
	if cur_seq_name is not None and not cur_seq:
	    raise Parse_error()
		
	# Output current sequence, if we're not at the start of the file
        if cur_seq_name is not None:
            yield (cur_seq_name, sequence_from_string(''.join(cur_seq)) )
            cur_seq = [ ]
            
        if not line: break
        
        cur_seq_name = line[1:].strip().split()[0]

def eland_iterator(filename):
    f = open(filename,'rU')
    while True:
        line = f.readline()
        if not line: break
        
        parts = line.split()
        
        if len(parts) < 2:
            raise Parse_error()
        
        yield (parts[0], sequence_from_string(parts[1]))

def sequence_file_iterator(filename):
    try:
        for result in fasta_iterator(filename):
            yield result
    except Parse_error:
        for result in eland_iterator(filename):
            yield result

def sequence_files_iterator(filenames):
    for filename in filenames:
        for result in sequence_file_iterator(filename):
            yield result
