#!/usr/bin/env python

"""
basic module to use as an example.
"""
import argparse


def add_two_numbers(a, b):
    """

    Args:
        a: int, first number to add
        b: int, second number to add
    Returns:
        c: int, addition between a and b
    """

    c = a + b
    return c


def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--firstnum",
        required=True,
        type=int,
        default=0
    )
    parser.add_argument(
        "--secondnum",
        required=True,
        type=int,
        default=0
    )
    args = parser.parse_args()

    first_num = args.firstnum
    second_num = args.secondnum

    added = add_two_numbers(a=first_num, b=second_num)
    print("{} + {} = {}".format(first_num, second_num, added))

if __name__ == "__main__":
    main()
