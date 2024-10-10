def compress_array(array):
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
    array_str_tuples = array_str[2:-2].split('), (')
    for tuple_str in array_str_tuples:
        x, y = tuple_str.split(', ')
        compressed_array.append((int(x), int(y)))
    return compressed_array