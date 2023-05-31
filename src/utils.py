# -*- coding: utf-8 -*-

import argparse

def get_args(*args) -> argparse.Namespace:
    """
    Wrapper for argparse.ArgumentParser. Receives a list of tuples,
    each tuple consisting of (<name_arg>, <default>) or (<name_arg>,),
    and passes the information contained within every tuple to an
    argument parser (instance of argparse.ArgumentParser). The
    arguments are, then, parsed (from sys.argv) and returned as an
    instance of argparse.Namespace.
    
    Parameters
    ----------
    args: list
        A list containing an arbitrary number of tuples
        (<name_arg>, <default>) / (<name_arg>,)
    """
    parser = argparse.ArgumentParser()
    for tup in args:
        name, *temp = tup
        default = temp[0] if temp else None
        parser.add_argument(name, default=default)
    args = parser.parse_args()
    return args
