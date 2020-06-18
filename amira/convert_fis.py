#! /usr/bin/env python3
# Convert output of ghrane's to the JSON format used by amira.
# The line "1 5 4 6 (245)" is transformed to "\"1_4_5_6\": 245".
# Commas and braces are added as appropriate to have a valid JSON output.

# Copyright 2019 Matteo Riondato <riondato@acm.org>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys

def main():
    with open(sys.argv[1]) as f:
        first = f.readline()
        print("{")
        sys.stdout.write("\t\"*\": {}".format(first.strip("()\n")))
        for line in f:
            o = line.index("(")
            its = [] # items
            for i in line[0:o-1].split(" "):
                its.append(int(i))
            its.sort()
            sys.stdout.write(",\n\t\"{}\": {}".format(
                "_".join(str(i) for i in its), line[o+1:-2]))
        sys.stdout.write("\n}\n")


if __name__ == '__main__':
    main()
