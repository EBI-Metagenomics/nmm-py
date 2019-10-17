def make_sure_bytes(p):
    try:
        p = p.encode()
    except AttributeError:
        pass
    return p
