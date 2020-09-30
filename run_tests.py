from testsuite.test_runner import TestRunner
import testsuite.utils as utils
import argparse
import os
import multiprocessing as mp
import itertools as itt


def get_parser():
    parser = argparse.ArgumentParser('tests the one-network-fits-all hypothesis')
    subparsers = parser.add_subparsers(dest='mode', required=True, help='run tests sequentially or in parallel')
    sequential = subparsers.add_parser('sequential', help='runs the tests sequentially')
    sequential.add_argument('--network', type=utils.GGINetworkSelector, choices=list(utils.GGINetworkSelector), required=True)
    sequential.add_argument('--generator', type=utils.NetworkGeneratorSelector, choices=list(utils.NetworkGeneratorSelector), required=True)
    sequential.add_argument('--verbose', action='store_true', help='print progress to stdout')
    parallel = subparsers.add_parser('parallel', help='runs the tests in parallel')
    return parser


def run_tests(ggi_network_selector, network_generator_selector, verbose=False):
    try:
        if verbose:
            print('loading data ...')
        test_runner = TestRunner()
        if verbose:
            print('running the tests ...')
        test_runner.run_all(ggi_network_selector, network_generator_selector, verbose)
        if verbose:
            print('saving the results')
        test_runner.save_results()
        return 0
    except:
        return 1


if __name__ == '__main__':
    os.chdir('testsuite')
    args = get_parser().parse_args()
    if args.mode == 'sequential':
        print('running tests sequentially ...')
        run_tests(args.network, args.generator, args.verbose)
    elif args.mode == 'parallel':
        print('running tests in parallel ...')
        networks = list(utils.GGINetworkSelector)
        generators = list(utils.NetworkGeneratorSelector)
        pool = mp.Pool(len(networks) * len(generators))
        exit_codes = pool.starmap(run_tests, list(itt.product(networks, generators)))
        print(f'exit codes: {exit_codes}')



