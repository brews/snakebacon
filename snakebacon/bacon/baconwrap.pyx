from libc.stdlib cimport malloc, free


def cook(str infile, str outfile, int ssize):
    cdef extern from "bacon.cpp":
        int notmain(int argc, char *argv[])
    cdef char **outgoing_argv
    cdef bytes binfile
    cdef bytes boutfile
    cdef bytes bssize
    cdef char* cinfile
    cdef char* coutfile
    cdef char* cssize
    cdef int argcount = 4
    outgoing_argv = <char**> malloc(argcount * sizeof(char*))
    try:
        binfile = infile.encode('utf8')
        cinfile = binfile

        boutfile = outfile.encode('utf8')
        coutfile = boutfile

        bssize = str(ssize).encode('utf8')
        cssize = bssize
        
        outgoing_argv[1] = cinfile
        outgoing_argv[2] = coutfile
        outgoing_argv[3] = cssize
        return notmain(argcount, outgoing_argv)
    finally:
        free(outgoing_argv)

# cook('MSB2K_20.bacon', 'out.bacon', 2000)