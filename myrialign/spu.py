
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
     This module contains a routine to compile C programs for the Cell SPU,
     and cache the result.
     
     It also contains SPU versions of routines found in other files in this
     package.

"""

import sys, os, os.path, sha, fcntl

import util

cache_dir = os.path.join(os.environ['HOME'],'.spucache')
lock_filename = os.path.join(cache_dir, 'lock') 

compile_command = 'spu-gcc -O5 -funroll-loops -o %(out_filename)s %(c_filename)s'

def get(code):
    if not os.path.isdir(cache_dir):
        os.mkdir(cache_dir)

    hasher = sha.new()
    hasher.update(compile_command)
    hasher.update(chr(0))
    hasher.update(code)
    filename = os.path.join(cache_dir, hasher.hexdigest())
    
    if os.path.isfile(filename):
        return filename
    
    lockfile = open(lock_filename,'ab')
    fcntl.lockf(lockfile, fcntl.LOCK_EX)
    try:
        if os.path.isfile(filename):
            return filename
        
        c_filename = os.path.join(cache_dir, 'thing.c')
        out_filename = os.path.join(cache_dir, 'thing')
        f = open(c_filename,'wb')
        f.write(code)
        f.close()
        
        #Hmmm
        util.show_status('Compiling helper')
        assert os.system(compile_command % locals()) == 0
        assert os.system('mv %(out_filename)s %(filename)s' % locals()) == 0
        return filename
    finally:
        fcntl.lockf(lockfile, fcntl.LOCK_UN)
            
    

matcher_defines = r"""
#define n_positions %(n_positions)d
#define n_errors %(n_errors)d
#define n_vecs %(n_vecs)d
#define indel_cost %(indel_cost)d
"""

matcher_body = r"""
#include <stdio.h>
#include <stdlib.h>
#include <spu_intrinsics.h>

#define BUFSIZE 8192

#define AT(i,j,k) ((i)*(n_positions*n_vecs)+(j)*n_vecs+(k))

static int position;
static vec_ullong2 match1[ n_errors*n_positions*n_vecs ], 
                   match2[ n_errors*n_positions*n_vecs ],
                   * __restrict__ matchin = match1,
                   * __restrict__ matchout = match2, 
                   nucmatches[ 5*n_positions*n_vecs ];

const vec_ullong2 ONES = { -1, -1 }, ZEROS = { 0, 0 };

static void error(char *error) {
    fprintf(stderr, "spu_match: %s\n", error);
    exit(1);
}

static void load(void *dest, size_t size, int n) {
    int result = fread(dest, size, n, stdin);
    if (result != n)
        error("Unexpected EOF");
}

static void output(int value) {
    fwrite(&value, sizeof(value), 1, stdout);
}

static void observe(unsigned char nuc) {
    int i, j, k, l, m;
    vec_ullong2 *temp, * __restrict__ this_nucmatches;
    
    this_nucmatches = nucmatches + (n_positions*n_vecs*nuc);
/*
// TODO: merge these
//  matchout[0,:] = nucmatches[:]
    for(j=0;j<n_positions;j++)
        for(k=0;k<n_vecs;k++)
            matchout[AT(0,j,k)] = this_nucmatches[j*n_vecs+k];

//  matchout[0,1:] &= matchin[0,:-1]
    for(j=1;j<n_positions;j++)
        for(k=0;k<n_vecs;k++)
            matchout[AT(0,j,k)] = 
               spu_and(matchout[AT(0,j,k)],
                       matchin [AT(0,j-1,k)]);
*/
    for(k=0;k<n_vecs;k++)
        matchout[k] = this_nucmatches[k];
    for(j=1;j<n_positions;j++)
        for(k=0;k<n_vecs;k++)
            matchout[AT(0,j,k)] = 
               spu_and(this_nucmatches[j*n_vecs+k],
                       matchin [AT(0,j-1,k)]);
        
//  for i in xrange(1,n_errors):
    for(i=1;i<n_errors;i++) {
//      matchout[i,i:] = nucmatches[i:]
        for(j=i;j<n_positions;j++)
            for(k=0;k<n_vecs;k++) {
                matchout[AT(i,j,k)] = this_nucmatches[j*n_vecs+k];

//      matchout[i,i+1:] &= matchin[i,i:-1]
//        for(j=i+1;j<n_positions;j++)
//            for(k=0;k<n_vecs;k++)
                if (j > i) {
                    matchout[AT(i,j,k)] = spu_and(
                        matchout[AT(i,j,k)],
                        matchin[AT(i,j-1,k)]
                    );
                }
                
                
//        for(j=i;j<n_positions;j++)
//            for(k=0;k<n_vecs;k++) {
//              #Mismatch
//              matchout[i,i:] |= matchin[i-1,i-1:-1]        
                matchout[AT(i,j,k)] = spu_or(
                    matchout[AT(i,j,k)],
                    matchin[AT(i-1,j-1,k)]
                );

//              # Deletion in read
//              matchout[i,i:] |= matchin[i-1,i:]
                if (i >= indel_cost)
                    matchout[AT(i,j,k)] = spu_or(
                        matchout[AT(i,j,k)],
                        matchin[AT(i-indel_cost,j,k)]
                    );
        
//              # Deletion in reference
//              matchout[i,i:] |= matchout[i-1,i-1:-1]
                if (i >= indel_cost)
                    matchout[AT(i,j,k)] = spu_or(
                        matchout[AT(i,j,k)],
                        matchout[AT(i-indel_cost,j-1,k)]
                    );
            }
    }
    
    // Any hits?
    for(k=0;k<n_vecs;k++) 
        for(l=0;l<2;l++) {
            unsigned long long a = spu_extract(matchout[AT(n_errors-1,n_positions-1,k)], l);
            if (__builtin_expect(a != 0, 0)) { // Hits are rare, don't expect them
                for(m=0;m<64;m++) {
                    unsigned long long mask = 1LLU<<(63-m); //Big endian 
                    if (a & mask) { 
                        //fprintf(stderr, "Hit found %d %d\n", position, k*128+l*64+m);
                        for(i=0;
                            i<n_errors &&
                            !(spu_extract(matchout[AT(i,n_positions-1,k)], l)&mask);
                            i++);
                            
                        //Superior prior match
                        //if (i &&
                        //    (spu_extract(matchin[AT(i-1,n_positions-1,k)], l)&mask) )
                        //   continue;

                        //Equivalent future match
                        //if (i &&
                        //    (spu_extract(matchout[AT(i-1,n_positions-2,k)], l)&mask) )
                        //   continue;
                           
                        //fprintf(stderr, "Hit found %d %d %d\n", position, k*128+l*64+m, i);
                        output(position);
                        output(k*128+l*64+m);
                        output(i);
                        fflush(stdout);
                    }
                }
            }
        }
    
    temp = matchin;
    matchin = matchout;
    matchout = temp;
    position += 1;
}

int main() {
    unsigned char buffer[BUFSIZE];
    int n_read, i,j,k;

    for(i=0;i<n_errors;i++)
        for(j=0;j<n_positions;j++)
            for(k=0;k<n_vecs;k++)
                match1[AT(i,j,k)] = 
                match2[AT(i,j,k)] = 
                    j < i ? ONES : ZEROS;
    
    //TODO: memory efficiency: don't actually store zeros for Ns
    load(nucmatches, sizeof(vec_ullong2), 5*n_positions*n_vecs);
    
    position = 0;
    
    while(n_read = fread(buffer, 1, BUFSIZE, stdin)) {
        int i;
        for(i=0;i<n_read;i++)
            observe(buffer[i]);
    } 
    
    return 0;
}
"""

def get_matcher(n_errors,n_positions,n_vecs,indel_cost):
    return get(matcher_defines % locals() + matcher_body)

#print get_matcher(4,5,6)

