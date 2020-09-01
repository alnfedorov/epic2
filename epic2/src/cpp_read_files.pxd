from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libc.stdint cimport uint32_t

ctypedef pair[string, char] key
ctypedef vector[uint32_t] intvec

ctypedef map[key, intvec] genome_map


cdef extern from "epic2/src/read_files.cpp":
    genome_map read_bed(const char *, uint32_t)
    genome_map read_bed_gz(const char *, uint32_t)
    genome_map read_bedpe(const char *, uint32_t)
    genome_map read_bedpe_gz(const char *, uint32_t)
    genome_map read_bam(const char *, uint32_t)
