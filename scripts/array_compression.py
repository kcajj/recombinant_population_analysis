import numpy as np


def compress_array(array):
    """
    Compress a numpy array by storing the number of consecutive occurances of each value.
    For example, the array [1, 2, 2, 3, 3, 3, 4, 4, 4, 4] would be compressed to [(1, 1), (2, 2), (3, 3), (4, 4)].
    :param array: The numpy array to compress
    :return: The compressed array
    """
    compressed_array = []
    previous = array[0]
    c = 0
    for i in range(len(array)):
        if array[i] == previous:
            c += 1
        else:
            compressed_array.append((c, int(previous)))
            c = 1
        previous = array[i]
    compressed_array.append((c, int(previous)))
    return compressed_array


def decompress_array(compressed_array):
    array = []
    for c, v in compressed_array:
        array += [int(v)] * c
    return array


def retrive_compressed_array_from_str(array_str):
    compressed_array = []
    array_str_tuples = array_str[2:-2].split("), (")
    for tuple_str in array_str_tuples:
        x, y = tuple_str.split(", ")
        compressed_array.append((int(x), int(y)))
    return compressed_array


def npz_extract(npz_file):
    npz = np.load(npz_file)
    lst = npz.files
    for item in lst:
        array = npz[item]
    return array
