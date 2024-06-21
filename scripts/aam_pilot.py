import argparse
import pathlib
import pandas as pd
import logging
import sys

def setup_logging(log_file):
    """Set up logging configuration."""
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger("AAMProcessing")
    return logger

def main(input_file, output_file, log_file, rsmi_column='reactions', mapper_types = ["rxn_mapper", "graphormer", "local_mapper"]):
    # Setup environment
    root_dir = pathlib.Path(__file__).resolve().parents[1]
    sys.path.append(str(root_dir))
    from SynTemp.SynAAM.atom_map_consensus import AAMConsensus
    from SynTemp.SynUtils.utils import load_database, save_database
    
    logger = setup_logging(log_file)
    logger.info(f"Starting processing of {input_file}")

    # Load data
    try:
        data = load_database(input_file)
    except Exception as e:
        logger.error(f"Failed to load database, trying to read CSV. Error: {e}")
        data = pd.read_csv(input_file).to_dict('records')

    # Mapping configurations
    mapper_types = ["rxn_mapper", "graphormer", "local_mapper"]
    rdt_jar_path = f"{root_dir}/Data/RDT_2.4.1.jar"
    working_dir = f"{root_dir}"

    # Process consensus mapping
    aam = AAMConsensus(data, mappers=mapper_types)
    logger.info(f"Starting consensus mapping for {input_file} using {mapper_types}")
    results = aam.batch_consensus(data, rsmi_column=rsmi_column, batch_size=len(data), job_timeout=None, safe_mode=False, rdt_jar_path=rdt_jar_path, working_dir=working_dir)
    logger.info(f"Completed consensus mapping for {input_file} using {mapper_types}")

    # Update data with results
    for i, entry in enumerate(data):
        for mapper in mapper_types:
            entry[mapper] = results[i].get(mapper, {})

    # Save results
    save_database(data, output_file)
    logger.info(f"Saved AAM results of {input_file} to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process AAM files and map atoms using different mapper types.")
    parser.add_argument("--input_file", type=str, required=True, help="Input file path for the AAM data.")
    parser.add_argument("--output_file", type=str, required=True, help="Output file path for the processed data.")
    parser.add_argument("--log_file", type=str, required=True, help="Log file path to record the processing progress.")
    parser.add_argument("--rsmi_column", type=str, default='reactions', help="Column name for RSMI data in the input file.")
    parser.add_argument("--mapper_types", type=list[str], default=["rxn_mapper", "graphormer", "local_mapper"], 
                        help="list of atom mapping techniques.")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.log_file, args.rsmi_column, args.mapper_types)
