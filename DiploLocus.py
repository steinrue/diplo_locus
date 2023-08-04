#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import argparse
import sys
# import os, subprocess


def main():
    parser = argparse.ArgumentParser(description="DiploLocus CLI")

    modes = parser.add_subparsers(title='DiploLocus Mode', dest='mode', metavar='',
                                  description='Subcommand to choose whether to compute likelihoods or simulate replicates of time-series samples.')
    # DL_lik.py
    from DiploLocus_likelihood import lik_args_parser

    lik_parser = modes.add_parser(
        name="likelihood",
        help="Computes loglikelihood for time-series samples of independent loci.\n"
    )
    lik_parser = lik_args_parser(lik_parser)

    # DL_sim.py
    from DiploLocus_simulate import sim_args_parser
    sim_parser = modes.add_parser(
        name="simulate",
        help="Simulate and plot time-series samples and frequency trajectories of"
             " independent replicates under diploid selection.\n"
    )

    sim_parser = sim_args_parser(sim_parser)

    # then parse stuff
    all_args = parser.parse_args(sys.argv[1:])

    if all_args.mode == 'likelihood':
        import DiploLocus_likelihood
        DiploLocus_likelihood.main(all_args)

    elif all_args.mode == 'simulate':
        import DiploLocus_simulate
        DiploLocus_simulate.main(all_args)

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
