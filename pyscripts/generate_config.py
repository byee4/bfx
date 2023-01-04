#!/usr/bin/env python3

import socket
import random
import argparse
import os
import json

def is_port_in_use(port, hostname):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex((hostname, port)) == 0

    
def return_open_port(hostname):
    start = 1025
    stop = 8887
    port_in_use = True
    while port_in_use: 
        rand_port = random.randrange(start, stop)
        port_in_use = is_port_in_use(rand_port, hostname)

    return rand_port


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--hostname",
        required=False,
        default=None
    )
    parser.add_argument(
        "--outfile",
        required=False,
        default="config.json"
    )
    
    # Process arguments
    args = parser.parse_args()
    hostname = args.hostname if args.hostname is not None else socket.gethostbyname(socket.gethostname())
    outfile = args.outfile
    
    # main func
    config_contents = {
        "shiny.port": return_open_port(hostname),
        "shiny.host": hostname
    }
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(config_contents, f, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()
