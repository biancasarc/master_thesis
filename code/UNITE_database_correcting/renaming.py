from pygbif import species
from Bio import SeqIO
from functools import lru_cache
import csv
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--unite_db", required=True)
parser.add_argument("--curated_db", required=True)

args = parser.parse_args()

UNITE_database = args.unite_db
curated_reference_database = args.curated_db

changed_records_txt = "changed_records.txt"
changed_sh_csv = "changed_sh.csv"
changes_summary_txt = "stats.txt"
output_file = "output.fasta"




# statistics counts all headers, keeps track of all the changes made in each taxonomic level
statistics={'total_records':0, 'all_levels_found_match':0, 's':0, 'g':0, 'f':0, 'o':0, 'c':0, 'p':0, 'k':0} 

# rank_map is a dictionary for the shortcuts used in UNITE
# This is neccesary to avoid ambiguity in the case some entries have other classification such as "SUPERCLASS", 
#  which would be abbreviated to s. 

rank_map = {
    "KINGDOM": "k",
    "PHYLUM": "p",
    "CLASS": "c",
    "ORDER": "o",
    "FAMILY": "f",
    "GENUS": "g",
    "SPECIES": "s"
} 


# Cached GBIF lookup — avoids redundant API calls for repeated species
@lru_cache(maxsize=None)
def gbif_lookup(s, k, p, c, o, f, g):
    return species.name_backbone(scientificName=s, kingdom=k, phylum=p,
                                  class_=c, order=o, family=f, genus=g,
                                  checklistKey="7ddf754f-d193-4cc9-b351-99906754a03b")


# taxonomy_corrector takes the header of a sequence as a string, and returns the corrected header as found in GBIF/COL
def taxonomy_corrector(tax_string):  
    statistics['total_records'] +=1

    #1. Split the string into smaller strings, one string for each taxonomic level
    parts = tax_string.split("|")           # split once and reuse
    levels = parts[1].split(';')
    ids = parts[::2]  #this list will store the ids and doi in positions 0 and 1
    #[::2] means that only every other element will be stored in the list, starting from 0

    #2. Create a dictionary, where the key is the taxonomic level and the value is the name
    #   Example: {'k':'Fungi'; 'p': 'Basidiomycota' .... }
    sh= {}
    for level in levels:
        if '__' in level:
            key, value = level.split('__')
            sh[key] = value
    

    #3. Look up the name of each taxonomic level on GBIF. Results are stored in a (quite complex) dictionary
    gbif_result = gbif_lookup(sh['s'], sh['k'], sh['p'], sh['c'], sh['o'], sh['f'], sh['g'])

    classification = gbif_result["classification"][1:]  #the last part is to avoid adding the DOMAIN, which COL generates
    for entry in classification:
        if entry.get("rank")=="SPECIES":    # If all taxonomy levels including species are found in GBIF db
            statistics['all_levels_found_match'] +=1    # Adding to statistics that all levels were found
            entry["name"]=entry["name"].replace(" ", "_")
    

    #gbif_dict is a clean dictionary version of with the GBIF/COL results. format is 
#   {'c': 'Agaricomycetes','f': 'Peniophoraceae', 'g': 'Peniophora','k': 'Fungi','o': 'Russulales','p': 'Basidiomycota','s': 'Peniophora_albobadia'}
    gbif_dict = {rank_map[entry["rank"]]: entry["name"] for entry in classification if entry["rank"] in rank_map}

    for rank in gbif_dict:
        if gbif_dict[rank] != sh[rank]:
            statistics[rank]+=1
    
    taxonomy_string = ";".join(
        f"{rank}__{gbif_dict.get(rank, sh.get(rank))}" for rank in rank_map.values()) #if rank exists in gbif_dict, return that. If not, return the one in the original sh
    
    return f"{ids[0]}|{taxonomy_string}|{ids[1]}"


    # Creating a dictionary with all unique genera/species in one SH.
     # Example: {'SH1281904.10FU': ['Thelephora_albomarginata'], 'SH1282297.10FU': ['Tomentella_sp'], 'SH1146664.10FU': ['Odontia_sp']}
    
def sh_species_dictionary(infile):
    sh_species={}
    for record in SeqIO.parse(infile, "fasta"):
        sh = record.id.split("|")[2]
        sp = (record.id.split("s__")[1]).split("|")[0]
        if sh in sh_species:
            if sp not in sh_species[sh]:
                sh_species[sh].append(sp)
        else:
            sh_species[sh]=[sp]
    return sh_species

# Function that returns a dictionary with the SH's that have conflicts.
# Conflicts represent either multiple species or multiple generas in one SH
# Output dictionary is formatted like:
#  {'multiple_genera': ['SH1284818.10FU', 'SH1323040.10FU'], 'multiple_sp': ['SH1327650.10FU']}

def sh_conflicts_tracker(sh_species):
    # Step 2: Identifying SHs with multiple species/generas (problematic SHs)
    problematic_sh={'multiple_genera':[], 'multiple_sp':[]}
    for sh, species_list in sh_species.items():
        genera=[s.split('_')[0] for s in species_list]
        species=[s.split('_')[1] for s in species_list if s.split('_')[1]!='sp']
        if (len(set(genera))) >1: 
            problematic_sh['multiple_genera'].append(sh)
        if (len(set(species))) >1:     
            problematic_sh['multiple_sp'].append(sh)
            
    return problematic_sh

curated_records = []
with open(curated_reference_database) as infile:
    for record in SeqIO.parse(infile, "fasta"):
        if "|refs" in record.id:
            curated_records.append(record.id.split("|")[1])

# Build species dictionaries during the first pass to avoid re-reading files
input_species_dictionary = {}
output_species_dictionary = {}

def update_sh_dict(d, record_id):
    sh = record_id.split("|")[2]
    sp = (record_id.split("s__")[1]).split("|")[0]
    if sh in d:
        if sp not in d[sh]:
            d[sh].append(sp)
    else:
        d[sh] = [sp]

with (open(UNITE_database) as infile, open(output_file, "w+") as outfile, 
 open(changed_records_txt, "w") as changed_records, open(changed_sh_csv, "w",newline='') as changed_sh):
    for record in SeqIO.parse(infile, "fasta"):
        update_sh_dict(input_species_dictionary, record.id)  # build input dict on first pass

        if record.id.split("|")[0] in curated_records:  # if record is manually curated, we keep the original header and do not change it. We also write it to the output file without changing it.
            SeqIO.write(record, outfile, "fasta")
            update_sh_dict(output_species_dictionary, record.id)  # curated records are unchanged

        else:
            newid = taxonomy_corrector(record.id)
            oldid=record.id
            if newid != record.id:
                print(f"{record.id} \n   changed to \n{newid} \n ---------- ",file=changed_records)

            record.id = newid
            record.description = newid
            SeqIO.write(record, outfile, "fasta")
            update_sh_dict(output_species_dictionary, newid)  # build output dict on first pass

    # No need to seek and re-parse files — dictionaries are already built above

    g_conflict_resolved=0
    s_conflict_resolved=0

    initial_conflicts = sh_conflicts_tracker(input_species_dictionary)
    output_conflicts = sh_conflicts_tracker(output_species_dictionary)

    g_conflict_resolved = set(initial_conflicts['multiple_genera'])-set(output_conflicts['multiple_genera'])
    s_conflict_resolved = set(initial_conflicts['multiple_sp'])-set(output_conflicts['multiple_sp'])
    s_and_g_conflict_resolved = g_conflict_resolved.intersection(s_conflict_resolved)
    if len(s_and_g_conflict_resolved)>0:
        g_conflict_resolved = g_conflict_resolved - s_and_g_conflict_resolved
        s_conflict_resolved = s_conflict_resolved - s_and_g_conflict_resolved
        
    data = []
    for sh in g_conflict_resolved:
        data.append({'SH': sh,
            'Old_Entries': ", ".join(input_species_dictionary[sh]),
            'New_Entries': ", ".join(output_species_dictionary[sh]),
            'Conflict_Type': "Different genera in SH"})

    for sh in s_conflict_resolved:
        data.append({'SH': sh,
            'Old_Entries': ", ".join(input_species_dictionary[sh]),
            'New_Entries': ", ".join(output_species_dictionary[sh]),
            'Conflict_Type': "Different species in SH"})

    for sh in s_and_g_conflict_resolved:
        data.append({'SH': sh,
            'Old_Entries': ", ".join(input_species_dictionary[sh]),
            'New_Entries': ", ".join(output_species_dictionary[sh]),
            'Conflict_Type': "Different species and genera in SH"})

    fieldnames = ['SH', 'Old_Entries', 'New_Entries', 'Conflict_Type']
    writer = csv.DictWriter(changed_sh, fieldnames = fieldnames)
    writer.writeheader()
    writer.writerows(data)



with open(changes_summary_txt, "w") as outfile:
    print(f"STATISTICS: Out of {statistics['total_records']} records, {statistics['all_levels_found_match']} have matches on GBIF on all taxonomical levels including species.", file=outfile)
    print(f"{statistics['s']} records had the species corrected by GBIF", file=outfile)
    print(f"{statistics['g']} records had the genus corrected by GBIF", file=outfile)
    print(f"{statistics['f']} records had the family corrected by GBIF", file=outfile)
    print(f"{statistics['o']} records had the order corrected by GBIF", file=outfile)
    print(f"{statistics['c']} records had the class corrected by GBIF", file=outfile)
    print(f"{statistics['p']} records had the phylum corrected by GBIF", file=outfile)
