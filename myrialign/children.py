
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

     This module provides a way to create child processes and communicate 
     with them.

"""

import sys, os, subprocess, fcntl, select, struct, cPickle

WRITERS = { }

class Error(Exception): pass

class Write_to_dead_child(Error): pass

class Child:
    def __init__(self, invokation):
        self.running = False
        self.subprocess = subprocess.Popen(invokation,
                                           stdin=subprocess.PIPE,
                                           stdout=subprocess.PIPE,
                                           close_fds=True)        
        fcntl.fcntl(self.subprocess.stdin, fcntl.F_SETFL, os.O_NONBLOCK)
        self.stdin = self.subprocess.stdin
        self.stdin_closed = False
        self.stdout = self.subprocess.stdout
        self.closed = False
        self.return_code = None

    def read(self, amount):
        return read(amount, self.stdout)
        
    def write(self, data):
        self._check_status()    
        write(data, self.stdin)

    def close_stdin(self):
        close(self.stdin)
        self.stdin_closed = True

    def send(self, item):
        self._check_status()
        send(item, self.stdin)
    
    def receive(self):
        return receive(self.stdout)
    
    def kill(self, signal=9):
        abort_write(self.stdin)
        os.kill(self.subprocess.pid, signal)
    
    def close(self):
        global N_CHILDREN
        self.stdout.close()
        if not self.stdin_closed:
            self.close_stdin()
        flush(self.stdin)
        if self.return_code is None:
            self.return_code = self.subprocess.wait()
        self.closed = True

    def _check_status(self):
        if self.return_code is None:
            self.return_code = self.subprocess.poll()
        if self.return_code is not None:
            raise Write_to_dead_child()

class Self_child(Child):
    def __init__(self, invokation=[sys.executable, sys.argv[0], 'child']):
        Child.__init__(self, invokation)


def wait(readers=[], timeout=None):
    """ Write any pending stuff. """
    reader_reverse_map = { }
    reader_filenos = [ ]
    for reader in readers:
        if isinstance(reader, Child):
            number = reader.stdout.fileno()
        elif isinstance(reader, int):
            number = reader
        else:
            number = reader.fileno()
        reader_filenos.append(number)
        reader_reverse_map[number] = reader
    
    read_ready = [ ]
    while reader_filenos or WRITERS:
        read_ready, write_ready, _ = select.select(reader_filenos, WRITERS.keys(), [], timeout)
        
        if not write_ready: break
        
        for writer in write_ready:
            queue = WRITERS[writer]
            if queue[0] is None:
                writer.close()
                del queue[0]
            else:
                pos, string, on_error = queue[0]
                try:
                    n = os.write(writer.fileno(),string[pos:pos+1024])
                    queue[0][0] += n
                    if queue[0][0] == len(string): 
                        del queue[0]
                except OSError, exception: #Most likely broken pipe
                    del queue[0]
                    on_error(exception)
            if not queue:
                del WRITERS[writer]  

    return [ reader_reverse_map[number] for number in read_ready ]


def do_stuff():
    """ Call this during long computations if write queues may not be
        empty. """
    wait(timeout=0)


def read(size, file=sys.stdin):
    """ Read data from a file.
    """
    
    data = [ ]
    remainder = size
    while remainder:
        wait([file])
        new_data = os.read(file.fileno(), remainder)
        if not new_data: break #EOF
        data.append(new_data)
        remainder -= len(new_data)
    
    return ''.join(data)


def default_write_error(exception):
    raise exception

    
def write(data, file=sys.stdout, on_error=default_write_error):
    """ Write data to a file.
    
        File must be flushed and set to non-blocking. 
        Call flush() before closing.
        """
    if file not in WRITERS:
        WRITERS[file] = [ ]
    WRITERS[file].append([0,data,on_error])
    do_stuff()


def close(file=sys.stdout):
    if file not in WRITERS:
        WRITERS[file] = [ ]
    WRITERS[file].append(None)
    do_stuff()


def abort_write(file=sys.stdout):
    if file in WRITERS:
        del WRITERS[file]


def flush(file):
    while file in WRITERS: 
        wait()


def flush_all():
    while WRITERS:
        wait()


def send(object, file=sys.stdout, on_error=default_write_error):
    pickled = cPickle.dumps(object, 2) # 2 == binary format
    write(struct.pack('<q', len(pickled)), file, on_error)
    write(pickled, file, on_error)


def receive(file=sys.stdin):
    length_pack = read(8, file)
    if len(length_pack) < 8: 
        raise EOFError()
    length = struct.unpack('<q', length_pack)[0]
    data = read(length, file)
    if len(data) < length: 
        raise EOFError()
    return cPickle.loads(data)



def test(argv):
    if sys.argv[1:] == ['child']:    
        send(receive() + ' back at you')
        
        try:
            receive()
        except EOFError:
            print >> sys.stderr, 'EOF raised, as expected'
        
        return 0
    else:
        kids = [ Self_child(), Self_child() ]
        
        kids[0].send('hello')
        kids[1].send('world')
        while kids:
            for child in wait(kids):
                print child.receive()
                kids.remove(child)
                child.close()
        
        return 0
        
if __name__ == '__main__':
    sys.exit( test(sys.argv) )

