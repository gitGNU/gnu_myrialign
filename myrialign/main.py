
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

def show_help():
    print
    print 'Usage: myr [command] ...'
    print
    print 'Commands:'
    print
    print '    align - align reads to a reference'
    print 
    print '    shred - generate simulated Illumina reads'
    print
    print 'Enter just "myr [command]" for help on that command'
    print 
    return 1

def main(argv):
    if len(argv) == 1:
        return show_help()
    
    command = argv[1]
    argv = argv[2:]
    
    if command == 'align':
        import align
        return align.main(argv)
    elif command == 'child':
        import align
        return align.child(argv)
    
    elif command == 'shred':
        import shred
        shred.main(argv)
    else:    
        return show_help()


