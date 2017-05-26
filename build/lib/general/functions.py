'''
Created on Sep 17, 2015

@author: brian
'''

import os


def insert_newlines(string, every=80):
    """ Makes a long string easier to look at, adding newlines to every
        x-amount of characters.
        
        Returns: a string with newlines.
    """
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)


def fill(my_string,desired_len,pos):
    """ extends the given string to be the specified length.
        
        Args:
            my_string: the original string
            desired_len: how long the new string will be
            pos: if set to 1, the spaces will fill from the left.
                if set to -1, the spaces will fill from the right.
        Returns:
            the new string with the new length.
    """
    spaces=[' ']*desired_len
    if pos==-1:
        return (my_string+''.join(spaces[len(my_string):]))
    elif pos==1:
        return ((''.join(spaces[len(my_string):])+my_string))


def which(program):
    """ Returns the program path if it is in the system path
        Copied from: 
        http://stackoverflow.com/questions/377017/test-if-executable-exists-in-pyscripts
        
        Also can use shutil.which() in python3.3
    """
    
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def split_files(infile, every, startwith=None):
    """ Splits a single file into (x) files ;)
    
        Usage: 
            infile = "/home/brian/Documents/Genomes/peregrine_falcon/all_chr.fa"
            split_files(infile, 1000, '>')
        
        Args:
            infile: single file
            every: number of entries or lines in each file
            startwith: what each entry begins with (helpful for fasta splitting)
    
    """
    i = 0 # counts number of entries
    j = 1 # incremental file counter
    filename, ext = os.path.splitext(infile)
    o = open(filename+"_{0}".format(j)+ext, 'w')
    with open(infile, 'r') as f:
        for line in f:
            if(line.startswith(startwith) or startwith==None):
                i = i + 1
            if(i < every):
                o.write(line)
            else:
                i = 0
                j = j + 1
                o.close()
                o = open(filename+"_{0}".format(j)+ext, 'w')
                o.write(line)
            
def main():
    pass
if __name__ == '__main__':
    main()