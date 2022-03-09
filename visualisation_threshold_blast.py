import os
import argparse
import glob
import plotly.io as pio
from textwrap import dedent

from common import visualisation
from common import utils
from common import makeTable

##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description="Threshold helper" )

general_option = parser.add_argument_group(title = "General input dataset options")
general_option.add_argument("-f",'--fasta_folder',
                            metavar="<folder>",
                            dest="fastaFolder",
                            help="Path to the folder with all the fasta of each groups: one fasta per group",
                            required=True,
                            )
general_option.add_argument("-i",'--input_fasta',
                            metavar="<file>",
                            dest="input_fasta",
                            help="The path to the fasta file that was used for the blast all vs all",
                            required=True,
                            )                               
general_option.add_argument("-b",'--blastfile',
                            metavar="<file>",
                            dest="blastFile",
                            help="Blast all vs all file",
                            required=True,
                            )
general_option.add_argument("-css",'--cssfile',
                            metavar="<file>",
                            dest="css",
                            help="CSS file",
                            required=True,
                            )                            
general_option.add_argument("-o",'--output',
                            default=None,
                            dest="output",
                            metavar='<OUTPUT_DIR>',
                            help="Using <OUTPUT_DIR> for output files (default: blastFile directory)",
                            )
general_option.add_argument("-filter",'--filter_length',
                            default=100,
                            dest="filter",
                            metavar='<SIZE>',
                            type=int,
                            help="Use a filter in the length of the alignment (default: 100)",
                            )
general_option.add_argument("-lcc",'--length_choice_cov',
                            default='mean',
                            dest="length_choice_cov",
                            metavar='<METHOD>',
                            help=dedent("""
                            Length used for percentage overlap calculation 
                                       between 2 sequences:
                                       'mean'=mean of the 2 lengths (default),
                                       'subject'=subject length, 'query'=query length,
                                       'shortest'=shortest length, 'longest'=longest length
                            """),
                            )
general_option.add_argument("-id",'--length_choice_id',
                            default='mean',
                            dest="length_choice_id",
                            metavar='<METHOD>',
                            help=dedent("""
                            Length used for percentage identity calculation 
                                       between 2 sequences:
                                       'mean'=mean of the 2 lengths (default),
                                       'subject'=subject length, 'query'=query length,
                                       'shortest'=shortest length, 'longest'=longest length
                                       'HSP'=HSP length
                            """),
                            )
##########################################################################

args = parser.parse_args()

if args.output:
    OUTPUT = args.output
else:
    OUTPUT = os.path.dirname(args.blastfile)

utils.create_folder(OUTPUT)

output_table = os.path.join(OUTPUT, "blast_summarized.tsv")
output_removed = os.path.join(OUTPUT, "blast_notinfasta.tsv")
output_html = os.path.join(OUTPUT, "report_threshold.html")

if os.path.isfile(output_removed):
    os.remove(output_removed)

##########################################################################

if args.length_choice_id in ["mean","subject", "query", "shortest", "longest", "HSP"]:
    option_id = args.length_choice_id
else:
    sys.exit('Error:: The option in length for percentage of identity should be: "mean","subject", "query", "shortest", "longest", "HSP"')

##########################################################################

if args.length_choice_cov in ["mean","subject", "query", "shortest", "longest"]:
    option_id = args.length_choice_id
else:
    sys.exit('Error:: The option in length for percentage of coverage should be: "mean","subject", "query", "shortest", "longest"')

##########################################################################

# List of fasta
list_fasta = glob.glob(os.path.join(args.fastaFolder, "*.fasta"))

# Get protein length
prot_dict = makeTable.get_protein_info(args.input_fasta)

# Infer families 
df_fam = makeTable.get_cluster_info(list_fasta)

# Create blast summarised file
makeTable.create_table_threshold(blast_out = args.blastFile, 
                       families = df_fam,
                       protein_dict = prot_dict,
                       output = output_table,
                       output_removed = output_removed,
                       length_treshold = args.filter)

# Choose the plotl template
pio.templates.default = "plotly"

# Get all the different dataframes
all_fam_file = output_table

# Create the 2D plot
tmp_plot2D = visualisation.scatter2D_plotly(all_fam_file)

# Create the 3D plot
tmp_plot3D = visualisation.scatter3D_plotly(all_fam_file)

# Merge the two plots and create the HTML report
visualisation.fig2html(tmp_plot2D, tmp_plot3D, output_html, args.css)

##########################################################################

