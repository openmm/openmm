from __future__ import print_function
import sys
from os.path import join, exists
from subprocess import call

def main():
    # pass any extra arguments to the first ctest invocation
    command = ['ctest', '--output-on-failure'] + sys.argv[1:]
    print(command)
    if call(command) == 0:
        return 0

    # load the log file containing the failed tests and reformat
    # it in a way so that ctest can take it as input to rerun
    # the failing tests
    log = join('Testing', 'Temporary', 'LastTestsFailed.log')
    assert exists(log)

    failed = []
    for line in open(log):
        failed.append(line.split(':')[0])

    with open('FailedTests.log', 'w') as f:
        print(','.join(x + ',' + x for x in failed), file=f)

    return call(['ctest', '--output-on-failure', '-I', 'FailedTests.log'])


if __name__ == '__main__':
    sys.exit(main())

