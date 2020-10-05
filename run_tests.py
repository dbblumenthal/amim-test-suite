from testsuite.test_runner import TestRunner
import testsuite.utils as utils
import argparse
import os
import multiprocessing as mp
import itertools as itt
import traceback



def get_parser():
    parser = argparse.ArgumentParser('tests the one-network-fits-all hypothesis')
    subparsers = parser.add_subparsers(dest='mode', required=True, help='run tests sequentially or in parallel')
    sequential = subparsers.add_parser('sequential', help='runs the tests sequentially')
    sequential.add_argument('--network', type=utils.GGINetworkSelector, choices=list(utils.GGINetworkSelector), required=True)
    sequential.add_argument('--generator', type=utils.NetworkGeneratorSelector, choices=list(utils.NetworkGeneratorSelector), required=True)
    sequential.add_argument('--method', type=utils.AlgorithmSelector, choices=list(utils.AlgorithmSelector), required=True)

    sequential.add_argument('--verbose', action='store_true', help='print progress to stdout')
    parallel = subparsers.add_parser('parallel', help='runs the tests in parallel')
    parallel.add_argument('--method', type=utils.AlgorithmSelector, choices=list(utils.AlgorithmSelector), required=True)
    return parser


def run_tests(ggi_network_selector, network_generator_selector, method_selector = "all", verbose=False):
    try:
        if verbose:
            print('loading data ...')
        test_runner = TestRunner()
        if verbose:
            print('running the tests ...')
            if method_selector != "all":
                print("testing only {0} method".format(method_selector))
            else:
                print("testing all methods")
        test_runner.run_all(ggi_network_selector, network_generator_selector, method_selector, verbose)
        if verbose:
            print('saving the results')
        test_runner.save_results()
        return 0
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    os.chdir('testsuite')
    args = get_parser().parse_args()
    if args.mode == 'sequential':
        print('running tests sequentially ...')
        if args.method != None:
            run_tests(args.network, args.generator, args.method, args.verbose)
        else:
            run_tests(args.network, args.generator, args.method, args.verbose)

    elif args.mode == 'parallel':
        print('running tests in parallel ...')
        networks = list(utils.GGINetworkSelector)
        generators = list(utils.NetworkGeneratorSelector)
        pool = mp.Pool(len(networks) * len(generators))
        if args.method != None:
            exit_codes = pool.starmap(run_tests, list(itt.product(networks, generators, [args.method])))
        else:
            exit_codes = pool.starmap(run_tests, list(itt.product(networks, generators, [args.method])))
        print(f'exit codes: {exit_codes}')



