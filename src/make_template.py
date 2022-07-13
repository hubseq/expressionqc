import json
MODULE = 'expressionqc'

mi_template_json = {'module_version': '00.00.00', 'program_name': 'expressionqc', 'program_subname': '', 'program_version': '0.11.5', 'compute': {'environment': 'aws', 'language': 'Python', 'language_version': '3.7', 'vcpus': 2, 'memory': 6000}, 'program_arguments': '-type STAR', 'program_input': [{'input_type': 'file', 'input_file_type': 'TXT', 'input_position': -1, 'input_prefix': '-i'}, {'input_type': 'file', 'input_file_type': 'CSV', 'input_position': -1, 'input_prefix': '-i'}, {'input_type': 'file', 'input_file_type': 'TAB', 'input_position': -1, 'input_prefix': '-i'}], 'program_output': [{'output_type': 'folder', 'output_file_type': '', 'output_position': 0, 'output_prefix': '-o'}], 'alternate_inputs': [], 'alternate_outputs': [], 'defaults': {"output_file": ""}}
with open(MODULE+'.template.json','w') as fout:
    json.dump(mi_template_json, fout)

io_dryrun_json = {'input': ['s3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny1.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny2.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny3.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny4.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny5.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny6.ReadsPerGene.out.tab'], 'output': ['s3://hubseq-data/test/rnaseq/run_test1/expressionqc/'],  'alternate_inputs': [], 'alternate_outputs': [], 'program_arguments': '-groups group1,group1,group1,group2,group2,group2', 'sample_id': MODULE+'_tiny_test', 'dryrun': ''}
io_json = {'input': ['s3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny1.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny2.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny3.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny4.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny5.ReadsPerGene.out.tab', 's3://hubseq-data/test/rnaseq/run_test1/rnastar/rnastar_test_tiny6.ReadsPerGene.out.tab'], 'output': ['s3://hubseq-data/test/rnaseq/run_test1/expressionqc/'],  'alternate_inputs': [], 'alternate_outputs': [], 'program_arguments': '-groups group1,group1,group1,group2,group2,group2', 'sample_id': MODULE+'_tiny_test'}

with open(MODULE+'.dryrun_test.io.json','w') as fout:
    json.dump(io_dryrun_json, fout)
with open(MODULE+'.test.io.json','w') as fout:
    json.dump(io_json, fout)

# job info test JSONs                                                                                                        
job_json = {"container_overrides": {"command": ["--module_name", MODULE, "--run_arguments", "s3://hubseq-data/modules/"+MODULE+"/job/"+MODULE+".test.job.json", "--working_dir", "/home/"]}, "jobqueue": "batch_scratch_queue", "jobname": "job_"+MODULE+"_test"}
with open(MODULE+'.test.job.json','w') as fout:
    json.dump(io_json, fout)
