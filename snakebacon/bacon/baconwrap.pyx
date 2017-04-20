from libc.stdlib cimport malloc, free


def cook(str infile, str outfile, int ssize):
    cdef extern from "bacon.c":
        int main(int argc, char *argv[])
    cdef char **outgoing_argv
    cdef bytes binfile
    cdef bytes boutfile
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
        outgoing_argv[1] = cinfile
        outgoing_argv[2] = coutfile
        outgoing_argv[3] = cssize
        return main(argcount, outgoing_argv)
    finally:
        free(outgoing_argv)
