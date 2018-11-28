"""
Run test suite through CTest, with some options set for the CI environment.

- Runs with a per-test and overall timeout which can be governed by the
  avialable time on the CI system.
- Reruns tests which fail (does not rerun tests which merely timeout).

"""
from __future__ import print_function
import sys
import os.path
import shutil
import time
from glob import glob
from os.path import join, exists
from subprocess import call
from argparse import ArgumentParser
from datetime import datetime, timedelta
from xml.etree import ElementTree


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--start-time",
        help="Time at which the overall CI job started (unix timestamp)",
        type=int,
        default=int(time.time()))
    parser.add_argument(
        "--job-duration",
        help="Overall time budget for the CI job (minutes). Default=30",
        default=30.,
        type=float)
    parser.add_argument(
        "--run-percent",
        help="Allocate this percent test execution time executing the main "
        "test suite. The remaining fraction will be used for re-running "
        "failing tests. Default=90",
        type=float,
        default=90.)
    parser.add_argument(
        '--timeout',
        help="Timeout for individual tests (seconds). Default=180",
        type=str,
        default='180')
    parser.add_argument(
        '--in-order',
        help='Run the tests in order',
        default=False,
        action='store_true')
    parser.add_argument(
        '--parallel',
        help='Number of processors to use',
        type=int,
        default=1)

    args, raw_args = parser.parse_known_args()

    status = execute_tests(args, raw_args)
    if status != 0:
        status = execute_failed_tests(args, raw_args)
    return status


def execute_tests(options, raw_options):
    start_time = datetime.fromtimestamp(options.start_time)
    stop_time = start_time + timedelta(minutes=options.job_duration)

    # timedelta for the amount of time from now until the CI job runs out
    remaining = stop_time - datetime.now()

    # tell CTest only to use some fraction of the remaining time for this
    # invocation
    stop_time = start_time + timedelta(
        seconds=(options.run_percent / 100.0) * remaining.seconds)

    if os.path.isdir('Testing'):
        shutil.rmtree('Testing')
    return call(['ctest',
                 '--output-on-failure',
                 '--parallel', str(options.parallel),
                 '-T', 'Test',
                 '--timeout', options.timeout,
                 '--stop-time', stop_time.strftime('%H:%M:%S')] + raw_options +
                 (['--schedule-random'] if options.in_order else []))


def execute_failed_tests(options, raw_options):
    matches = glob('Testing/*/Test.xml')
    assert len(matches) == 1
    root = ElementTree.parse(matches[0])
    tests = root.findall('.//Testing/Test')

    def failed_without_timeout(test_node):
        if test_node.get('Status') == 'failed':
            return test_node.find(
                'Results/NamedMeasurement[@name="Exit Code"]/Value').text != 'Timeout'

    failed_tests = [t.find('Name').text
                    for t in tests if failed_without_timeout(t)]
    print('*'*30)
    print('Rerunning failing tests...')
    print('*'*30)

    start_time = datetime.fromtimestamp(options.start_time)
    stop_time = start_time + timedelta(minutes=options.job_duration)
    return call(['ctest'] + raw_options + [
                 '--output-on-failure',
                 '--parallel', str(options.parallel),
                 '-R', '|'.join(failed_tests),
                 '--timeout', options.timeout,
                 '--stop-time', stop_time.strftime('%H:%M:%S')] +
                 (['--schedule-random'] if options.in_order else []))


if __name__ == '__main__':
    sys.exit(main())
