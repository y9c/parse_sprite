cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import re

cimport cedlib



from libc.stdlib cimport malloc


cdef char *basemap = [ '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',   '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a' ]

def reverse_complement(str seq):
    cdef int seq_len = len(seq)
    cdef char *seq_dest = <char *>malloc(seq_len + 1)
    seq_dest[seq_len] = '\0'

    cdef bytes py_bytes = seq.encode('UTF-8')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('UTF-8')


class NeedsAlphabetMapping(Exception):
    pass


def _map_ascii_string(s):
    """ Helper func to handle the bytes mapping in the case of ASCII strings."""
    if isinstance(s, bytes):
        return s
    elif isinstance(s, str):
        s_bytes = s.encode('utf-8')
        if len(s_bytes) == len(s):
            return s_bytes
    raise NeedsAlphabetMapping()


def _map_to_bytes(query, target, additional_equalities):
    """ Map hashable input values to single byte values.

    Example:
    In: query={12, 'ä'}, target='ööö', additional_equalities={'ä': 'ö'}
    Out: b'\x00\x01', b'\x02\x02\x02', additional_equalities={b'\x01': b'\x02'}
    """
    cdef bytes query_bytes
    cdef bytes target_bytes
    try:
        query_bytes = _map_ascii_string(query)
        target_bytes = _map_ascii_string(target)
    except NeedsAlphabetMapping:
        # Map elements of alphabet to chars from 0 up to 255, so that Edlib can work with them,
        # since C++ Edlib needs chars.
        alphabet = set(query).union(set(target))
        if len(alphabet) > 256:
            raise ValueError(
                "query and target combined have more than 256 unique values, "
                "this is not supported.")
        alphabet_to_byte_mapping = {
            c: idx.to_bytes(1, byteorder='big') for idx, c in enumerate(alphabet)
        }
        map_seq = lambda seq: b''.join(alphabet_to_byte_mapping[c] for c in seq)
        query_bytes = map_seq(query)
        target_bytes = map_seq(target)
        if additional_equalities is not None:
            additional_equalities = [
                (alphabet_to_byte_mapping[a].decode('utf-8'), alphabet_to_byte_mapping[b].decode('utf-8'))
                for a, b in additional_equalities
                if a in alphabet_to_byte_mapping and b in alphabet_to_byte_mapping]
    return query_bytes, target_bytes, additional_equalities


def align(query, target, mode="NW", task="distance", k=-1, additionalEqualities=None):
    """ Align query with target using edit distance.
    @param {str or bytes or iterable of hashable objects} query, combined with target must have no more
           than 256 unique values
    @param {str or bytes or iterable of hashable objects} target, combined with query must have no more
           than 256 unique values
    @param {string} mode  Optional. Alignment method do be used. Possible values are:
            - 'NW' for global (default)
            - 'HW' for infix
            - 'SHW' for prefix.
    @param {string} task  Optional. Tells edlib what to calculate. The less there is to calculate,
            the faster it is. Possible value are (from fastest to slowest):
            - 'distance' - find edit distance and end locations in target. Default.
            - 'locations' - find edit distance, end locations and start locations.
            - 'path' - find edit distance, start and end locations and alignment path.
    @param {int} k  Optional. Max edit distance to search for - the lower this value,
            the faster is calculation. Set to -1 (default) to have no limit on edit distance.
    @param {list} additionalEqualities  Optional.
            List of pairs of characters or hashable objects, where each pair defines two values as equal.
            This way you can extend edlib's definition of equality (which is that each character is equal only
            to itself).
            This can be useful e.g. when you want edlib to be case insensitive, or if you want certain
            characters to act as a wildcards.
            Set to None (default) if you do not want to extend edlib's default equality definition.
    @return Dictionary with following fields:
            {int} editDistance  Integer, -1 if it is larger than k.
            {int} alphabetLength Integer, length of unique characters in 'query' and 'target'
            {[(int, int)]} locations  List of locations, in format [(start, end)].
            {string} cigar  Cigar is a standard format for alignment path.
                Here we are using extended cigar format, which uses following symbols:
                Match: '=', Insertion to target: 'I', Deletion from target: 'D', Mismatch: 'X'.
                e.g. cigar of "5=1X1=1I" means "5 matches, 1 mismatch, 1 match, 1 insertion (to target)".
    """
    # Transform python sequences of hashables into c strings.
    cdef bytes query_bytes
    cdef bytes target_bytes
    query_bytes, target_bytes, additionalEqualities = _map_to_bytes(
            query, target, additionalEqualities)
    cdef char* cquery = query_bytes;
    cdef char* ctarget = target_bytes;

    # Build an edlib config object based on given parameters.
    cconfig = cedlib.edlibDefaultAlignConfig()

    if k is not None: cconfig.k = k

    if mode == 'NW': cconfig.mode = cedlib.EDLIB_MODE_NW
    if mode == 'HW': cconfig.mode = cedlib.EDLIB_MODE_HW
    if mode == 'SHW': cconfig.mode = cedlib.EDLIB_MODE_SHW

    if task == 'distance': cconfig.task = cedlib.EDLIB_TASK_DISTANCE
    if task == 'locations': cconfig.task = cedlib.EDLIB_TASK_LOC
    if task == 'path': cconfig.task = cedlib.EDLIB_TASK_PATH

    cdef bytes tmp_bytes;
    cdef char* tmp_cstring;
    cdef cedlib.EdlibEqualityPair* c_additionalEqualities = NULL
    if additionalEqualities is None:
        cconfig.additionalEqualities = NULL
        cconfig.additionalEqualitiesLength = 0
    else:
        c_additionalEqualities = <cedlib.EdlibEqualityPair*> PyMem_Malloc(len(additionalEqualities)
                                                                          * cython.sizeof(cedlib.EdlibEqualityPair))
        for i in range(len(additionalEqualities)):
            c_additionalEqualities[i].first = bytearray(additionalEqualities[i][0].encode('utf-8'))[0]
            c_additionalEqualities[i].second = bytearray(additionalEqualities[i][1].encode('utf-8'))[0]
        cconfig.additionalEqualities = c_additionalEqualities
        cconfig.additionalEqualitiesLength = len(additionalEqualities)

    # Run alignment -- need to get len before disabling the GIL
    query_len = len(query)
    target_len = len(target)
    with nogil:
        cresult = cedlib.edlibAlign(cquery, query_len, ctarget, target_len, cconfig)
    if c_additionalEqualities != NULL: PyMem_Free(c_additionalEqualities)

    if cresult.status == 1:
        raise Exception("There was an error.")

    # Build python dictionary with results from result object that edlib returned.
    locations = []
    if cresult.numLocations >= 0:
        for i in range(cresult.numLocations):
            locations.append((cresult.startLocations[i] if cresult.startLocations else None,
                              cresult.endLocations[i] if cresult.endLocations else None))
    cigar = None
    if cresult.alignment:
        ccigar = cedlib.edlibAlignmentToCigar(cresult.alignment, cresult.alignmentLength,
                                              cedlib.EDLIB_CIGAR_EXTENDED)
        cigar = <bytes> ccigar
        cigar = cigar.decode('UTF-8')
    result = {
        'editDistance': cresult.editDistance,
        'alphabetLength': cresult.alphabetLength,
        'locations': locations,
        'cigar': cigar
    }
    cedlib.edlibFreeAlignResult(cresult)

    return result


def get_min_edit_distance(query, refs, min_dist=2):
    """Get the minimum edit distance between two sequences."""
    for n, s in refs.items():
        dist = align(query, s, mode="NW")["editDistance"]
        if dist <= min_dist:
            return n
    return "NA"



def overlap_pairend(seq1, seq2, min_overlap=3, left_shift=9, right_shift=22):
    """seq2 has been reverse complemented."""
    for start in range(left_shift, len(seq1) - min_overlap):
        l = min(len(seq1) - start, len(seq2))
        dist = align(seq1[start : start + len(seq2)], seq2[:l], mode="NW")[
            "editDistance"
        ]
        if dist / l < 0.1:
            return max(start + right_shift, 0)
    return None
