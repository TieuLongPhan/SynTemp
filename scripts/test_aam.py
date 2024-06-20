import pathlib
from SynTemp.SynAAM.atom_map_consensus import AAMConsensus
root_dir = pathlib.Path(__file__).parents[0]


mapper_types = ["rxn_mapper", "graphormer", "local_mapper", "rdt"]

rdt_jar_path = f"{root_dir}/Data/RDT_2.4.1.jar"
working_dir = f"./"

data = [{'id':1, 'reactions': "CCOP(=O)(CC(=O)O)OCC.CC=O>>CCOP(=O)(O)OCC"},
        {'id':2, 'reactions': "CCOP(=O)(CC(=O)O)OCC.CC=O.[HH]>>CCOP(=O)(O)OCC.CCCC(=O)O"}]

aam = AAMConsensus(data, mappers=mapper_types)

results = aam.batch_consensus(data, rsmi_column='reactions', batch_size=len(data), job_timeout=None,
                              safe_mode=False, rdt_jar_path=rdt_jar_path, working_dir=working_dir)

print(results)