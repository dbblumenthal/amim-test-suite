from testsuite.test_runner import TestRunner
import testsuite.utils as utils
import argparse

def get_parser():
    parser = argparse.ArgumentParser('tests the one-network-fits-all hypothesis')
    parser.add_argument('--network', type=utils.GGINetworkSelector, choices=list(utils.GGINetworkSelector), required=True)
    parser.add_argument('--generator', type=utils.NetworkGeneratorSelector, choices=list(utils.NetworkGeneratorSelector), required=True)
    parser.add_argument('--verbose', action='store_true', help='print progress to stdout')
    return parser


if __name__ == '__main__':
    args = get_parser().parse_args()
    test_runner = TestRunner()
    test_runner.run_all(args.network, args.num_randomizations, args.verbose)
    test_runner.save_results()