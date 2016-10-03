import argparse
import genetic_analyzer as ga
import sys


def main():
    # Parse the command line inputs
    parser = argparse.ArgumentParser(description="profile_normalization")
    parser.add_argument("-settings", dest="settings", required=True, help="settings.txt", metavar="string")
    # parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
    args = parser.parse_args()
    # Run the command
    # samples = ga.load_settings(args.samples)
    settings = ga.load_settings(args.settings)
    sample = ""
    for i in settings:
        if i != None:
            sample = i
            continue
    promoter_fn = ga.promoter_reu_filename(settings, sample)
    terminator_fn = ga.terminator_reu_filename(settings, sample)

    # Merge promoter and terminator expected expression data with observed pipeline data
    ga.merge_reus(promoter_fn, ga.combined_promoter_profile_perf_filename(settings))
    ga.merge_reus(terminator_fn, ga.combined_terminator_profile_perf_filename(settings))

    # Set new outfile name for promoter images
    part_perf_reu_name = ga.combined_promoter_profile_perf_filename(settings)[:-3] + 'reu.txt'
    print("Making promotor comparison graphs", part_perf_reu_name)
    status1 = ga.comparison_graphs(settings, part_perf_reu_name)

    # Set new outfile name for terminator images
    part_perf_reu_name = ga.combined_terminator_profile_perf_filename(settings)[:-3] + 'reu.txt'
    print("Making terminator comparison graphs", part_perf_reu_name)

    status2 = ga.comparison_graphs(settings, part_perf_reu_name)

    return status1 + status2


if __name__ == "__main__":
    status = main()
    sys.exit(status)
