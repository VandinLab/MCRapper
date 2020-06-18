#! /usr/bin/env python3
import re
import sys

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Error: wrong number of arguments\nUsage: {} "
                "FILE\n".format(sys.argv[0]))
        sys.exit(1)
    with open(sys.argv[1], 'rt') as f:
        reg = re.compile("^\s+\"(mine|get_rho2)\":\s[0-9]+$")
        for line in f:
            res = reg.search(line)
            if res:
                sys.stdout.write("{},\n".format(
                    res.group(0)))
            else:
                sys.stdout.write(line)

if __name__ == "__main__":
    main()
