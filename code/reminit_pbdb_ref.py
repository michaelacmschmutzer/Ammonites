# Remove all the initials from reflist

import re

with open('../data/pbdb_reflist.txt', 'r') as f:
    lines = f.readlines()

new_lines = []
for l in lines:
    # Remove all combinations of an uppercase letter followed by a dot
    lin = re.sub(r'[A-Z]\.', '', l)
    # Clear spaces in front of first author last name
    lfs = re.sub(r'\s+', '', lin, count = 1)
    # Clear all multiple spaces
    las = re.sub(r'\s{2,}', ' ', lfs)
    new_lines.append(las)

with open('../data/pbdb_reflist_noinit.txt', 'w') as f:
    for line in new_lines:
        f.write(line)
