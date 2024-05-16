from pathlib import Path

# CONSTANTS

# paths
PROJECTS_PATH = Path("projects/")
BATCH_SCRIPTS_PATH = Path("batch_scripts/")
SCRIPTS_PATH = Path("scripts/")

# configs
DEFAULT_MODE = 'default'
DEFAULT_DATABASE = 'test'

# Scripts to be run in the default mode
SCRIPTS = {
    1: '01_search.sh',
    2: '02_mk_search_df.py',
    3: '03_plot_search.py',
    4: '04_filter.py',
    5: '05_cluster.sh',
    6: '06_mk_genome_consistent.py',
    7: '07_cluster_dict.py',
    8: '08_plot_filtering_clustering.py',
    9: '09_align.sh',
    10: '10_get_translations.py',
    11: '11_get_gen_context.py',
    12: '12_trim.sh',
    13: '13_trim_report.py',
    14: '14_domain_search.sh',
    15: '15_prottest.sh',
    16: '16_domain_data.py',
    17: '17_tree.sh'
}

# Scripts to be run in the cluster mode
BATCH_SCRIPTS = {
    1: '01_blastp.batch',
    2: '02_mk_blastp_df.batch',
    3: '03_plot_blastp.batch',
    4: '04_filter_blastp.batch',
    5: '05_cluster.batch',
    6: '06_mk_genome_consistent.batch',
    7: '07_cluster_dict.batch',
    8: '08_plot_filtering_clustering.batch',
    9: '09_align.batch',
    10: '10_get_translations.batch',
    11: '11_get_gen_context.batch',
    12: '12_trim.batch',
    13: '13_trim_report.batch',
    14: '14_hmmscan.batch',
    15: '15_prottest.batch',
    16: '16_get_domains.batch',
    17: '17_raxml.batch'
}

# Prompt listing available commands
PROMPT = '''
Prompt:
"n [project_name] [query_path]": create new project [project_name] with query at [query_path]
"s [project] [step]": start [project] [step]
"a": start all ready steps in all projects
"d [project]": delete [project]
"q": quit
'''

# This is written into the project status file `status.txt` when project is created
STATUS_FILE = '''\t{}
1\tready
2\t-
3\t-
4\t-
5\t-
6\t-
7\t-
8\t-
9\t-
10\t-
11\t-
12\t-
13\t-
14\t-
15\t-
16\t-
17\t-
'''

# Steps of the analysis. Used to create status dataframe
STEPS = {'Step': ['1. Search',
                  '2. Hit data',
                  '3. Search report',
                  '4. Filter',
                  '5. Cluster',
                  '6. Genome consistence',
                  '7. Cluster data',
                  '8. Filtering and clustering report',
                  '9. Align',
                  '10. Full translations',
                  '11. Genome context',
                  '12. Trim',
                  '13. Trim report',
                  '14. Domain search',
                  '15. Prottest',
                  '16. Domain data',
                  '17. Tree']
         }

# Dependencies of the steps on each other
UNLOCKS = {
    1: [2],
    2: [3, 4],
    3: [],
    4: [5],
    5: [6, 7],
    6: [9, 10, 11],
    7: [8],
    8: [],
    9: [12],
    10: [13],
    11: [],
    12: [13, 14],
    13: [],
    14: [15],
    15: [16],
    16: [],
    17: []
}

