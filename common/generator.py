##########################################################################

def _make_gen(reader):
    """Generator that yields bytes from a reader .

    Args:
        reader (open file): The opened file that you want to read

    Yields:
        string??: a chuck part of the file
    """

    while True:
        b = reader(2 ** 16)
        if not b: break
        yield b

##########################################################################

def buf_count_newlines_gen(fname):
    """Generate the number of newlines in the given file .

    Args:
        fname (string): name of the file to read

    Returns:
        int: number of line in the file
    """

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count


##########################################################################

def buf_count_prot_gen(fname):
    """Generate the number of '>' in the given file .

    Args:
        fname (string): name of the file to read

    Returns:
        int: number of line in the file
    """

    with open(fname, "rb") as f:
        count = sum(buf.count(b">") for buf in _make_gen(f.raw.read))
    return count

##########################################################################