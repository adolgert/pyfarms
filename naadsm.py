from default_parser import DefaultArgumentParser



def load_naadsm_scenario(scenario_file, herd_file):
    


if __name__ == "__main__":
    parser=DefaultArgumentParser(description="Load NAADSM XML Files")
    parser.add_function("load", "Load XML files")
    parser.add_argument("--scenario", dest=scenario_file, type=string,
        action="store", default=None, help="Scenario file")
    parser.add_argument("--herd", dest=herd_file, type=string,
        action="store", default=None, help="Herd file")

    args=parser.parse_args()
    if args.load:
        load_naadsm_scenario(args.scenario_file, args.herd_file)
