
import sha, os, shutil, fcntl, sys

# TODO: spu should use this mechanism

cache_dir = os.path.join(os.environ['HOME'],'.myrcache')

def file_signature(filename):
    filename = os.path.abspath(filename)
    return (filename, os.path.getmtime(filename))

def get(ident, callback):
    """
        ident uniquely identifies this job
	should be a string or tuple of idents
        
    """
    # Bend over backwards to avoid race condition
    try:
        os.mkdir(cache_dir)
    except OSError:
        if not os.path.isdir(cache_dir):
	    raise
	
    hasher = sha.new()
    hasher.update(repr(ident))
    root = os.path.join(cache_dir, hasher.hexdigest())
    
    print >> sys.stderr, 'Cached job'
    print >> sys.stderr, ident
    print >> sys.stderr, root
    
    working_dir = root + '-working'
    lock_filename = root + '-lock'
    result_dir = root + '-result'
    
    lockfile = open(lock_filename,'ab')
    fcntl.lockf(lockfile, fcntl.LOCK_EX)
    try:
        if os.path.exists(result_dir):
            return result_dir

        if os.path.exists(working_dir):
	    shutil.rmtree(working_dir)
	os.mkdir(working_dir)
	    
	try:
	    callback(working_dir)
	except:
	    shutil.rmtree(working_dir)
	    raise
	
	os.rename(working_dir, result_dir)
	
        return result_dir
    finally:
        fcntl.lockf(lockfile, fcntl.LOCK_UN)

    return result_dir


def cached(function):
    """ Decorate a function with cachey goodness 
    
        Function will be passed a directory to place results in
        After decoration, function will return results directory
    """
    name = function.__name__
    module_name = function.__module__
    
    def inner(*args):
        return get((name,module_name,args), lambda working_dir: function(working_dir, *args))
    return inner
        



