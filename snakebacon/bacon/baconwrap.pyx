from libc.stdlib cimport malloc, free


def cook(args):
    cdef extern from "bacon.c":
            int main(int argc, char *argv[])
    cdef char **outgoing_argv
        # we put this in try/finally to assure we don't leak memory
    try:
        # first, allocate the memory.
        outgoing_argv = <char**> malloc(len(args) * sizeof(char*))
        # then we iterate through args to put the strings in our
        # char array, using Python C api to get char* from str
        for i in range(len(args)):
            outgoing_argv[i] = args[i]
            # outgoing_argv[i] = PyUnicode_AsUTF8(args[i])
        main(len(args), outgoing_argv)
    finally:
         # lastly, we want to make sure we free any
         # allocated memory.
        free(outgoing_argv)

    # # args = [bytes(x) for x in args]
    # cdef const char** c_argv = <const char**> malloc(len(args) * sizeof(char*))
    # # c_argv = <char**>malloc(sizeof(char*) * len(args))
    # if c_argv is NULL:
    #     raise MemoryError()
    # try:
    #     for idx, s in enumerate(args):
    #        c_argv[idx] = s
    #     # test(len(args), c_argv)
    # finally:
    #     free(c_argv)
    # main(len(args), c_argv)
    # main(2, "eggs")
