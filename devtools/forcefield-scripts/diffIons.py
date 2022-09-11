import sys
import re
import subprocess

file1, file2 = sys.argv[1], sys.argv[2]

with open(file1) as fh:
    text1 = fh.read()

with open(file2) as fh:
    text2 = fh.read()

m = re.search(r'Type class="([^\-]*)', text1)
if m:
    prefix1 = m.group(1)
m = re.search(r'Type class="([^\-]*)', text2)
if m:
    prefix2 = m.group(1)

text2a = text2.replace(prefix2, prefix1)
text2b = re.sub(r'<!--.*?-->', '', text2a)

with open('new.xml', 'w') as fh:
    fh.write(text2b)

subprocess.call(['diff', '-uw', file1, 'new.xml'])
