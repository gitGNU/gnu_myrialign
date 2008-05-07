
import sys

def show_status(message):
    sys.stderr.write('\r\x1b[K\r')
    sys.stderr.write(message)
    sys.stderr.flush()

def show_message(message):
    show_status(message+'\n')
