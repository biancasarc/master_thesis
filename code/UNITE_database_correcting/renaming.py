"""
-- CODE USAGE --
python renaming.py --unite_db input_UNITE_db.fasta --curated_db curated_UNITE_db.fasta

"""

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
modified_SH_csv = "modified_SH.csv"
summary_of_changes_txt = "stats.txt"
output_file = "output.fasta"



# statistics keeps track of all the changes made in each taxonomic level
statistics={'total_records':0, 'all_levels_found_match':0, 's':0, 'g':0, 'f':0, 'o':0, 'c':0, 'p':0, 'k':0} 

# rank_map is a dictionary for the shortcuts used in UNITE
# neccesary to avoid ambiguity if classifications like "SUPERCLASS" appear 

rank_map = {
    "KINGDOM": "k",
    "PHYLUM": "p",
    "CLASS": "c",
    "ORDER": "o",
    "FAMILY": "f",
    "GENUS": "g",
    "SPECIES": "s"
} 


# cached GBIF lookup — avoids redundant API calls for repeated species
@lru_cache(maxsize=None)
def gbif_lookup(s, k, p, c, o, f, g):
    return species.name_backbone(scientificName=s, kingdom=k, phylum=p,
                                  class_=c, order=o, family=f, genus=g,
                                  checklistKey="7ddf754f-d193-4cc9-b351-99906754a03b")



# taxonomy_corrector takes the header of a sequence and returns the corrected header as found in GBIF/COL
def taxonomy_corrector(tax_string):  
    statistics['total_records'] +=1

    # 1. Split the string - one string for each taxonomic level
    parts = tax_string.split("|")          
    levels = parts[1].split(';')
    id_and_doi = parts[::2]


    # 2. Create a dictionary, where the key is the taxonomic level and the value is the name
    #   Example: {'k':'Fungi'; 'p': 'Basidiomycota' .... }
    original_taxonomy= {}
    for level in levels:
        if '__' in level:
            key, value = level.split('__')
            original_taxonomy[key] = value
    

    # 3. Look up the name of each level on GBIF
    gbif_taxonomy = gbif_lookup(original_taxonomy['s'], original_taxonomy['k'], original_taxonomy['p'], 
                              original_taxonomy['c'], original_taxonomy['o'], original_taxonomy['f'], 
                              original_taxonomy['g'])

    classification = gbif_taxonomy["classification"][1:]
    
    # 4. Process the gbif output into the UNITE format 

    for entry in classification:
        if entry.get("rank")=="SPECIES":    
            statistics['all_levels_found_match'] +=1    # Adding to statistics that all levels were found
            entry["name"]=entry["name"].replace(" ", "_")
    

    #gbif_dict is a clean dictionary version of with the GBIF/COL results. format is 
    #   {'c': 'Agaricomycetes','f': 'Peniophoraceae', 'g': 'Peniophora','k': 'Fungi','o': 'Russulales','p': 'Basidiomycota','s': 'Peniophora_albobadia'}
    gbif_dict = {rank_map[entry["rank"]]: entry["name"] for entry in classification if entry["rank"] in rank_map}

    for rank in gbif_dict:
        if gbif_dict[rank] != original_taxonomy[rank]:
            statistics[rank]+=1
    
    output_string = ";".join(
        f"{rank}__{gbif_dict.get(rank, original_taxonomy.get(rank))}" for rank in rank_map.values()) #if rank exists in gbif_dict, return that. If not, return the one in the original sh
    
    return f"{id_and_doi[0]}|{output_string}|{id_and_doi[1]}"




# Output of sh_conflicts_tracker is formatted like:
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


# getting the SH's that have been manually curated 
curated_records = []
with open(curated_reference_database) as infile:
    for record in SeqIO.parse(infile, "fasta"):
        if "|refs" in record.id:
            curated_records.append(record.id.split("|")[1])



# add_species_to_SH updates a dictionary that tracks which species belong to each SH 
def add_species_to_SH(dictionary, record_id):
    sh = record_id.split("|")[2]
    sp = (record_id.split("s__")[1]).split("|")[0]
    if sh in dictionary:
        if sp not in dictionary[sh]:
            dictionary[sh].append(sp)
    else:
        dictionary[sh] = [sp]



input_species_dictionary = {}
output_species_dictionary = {}



with (open(UNITE_database) as infile, open(output_file, "w+") as outfile, 
 open(changed_records_txt, "w") as changed_records, 
 open(modified_SH_csv, "w",newline='') as modified_SH):
    

    for record in SeqIO.parse(infile, "fasta"):
        add_species_to_SH(input_species_dictionary, record.id)  #each species is added to the dictionary where SH is key

        if record.id.split("|")[0] in curated_records:  # if record is manually curated, it is written in the output.fasta as it is.
            SeqIO.write(record, outfile, "fasta")
            add_species_to_SH(output_species_dictionary, record.id)  

        else:
            newid = taxonomy_corrector(record.id)
            oldid=record.id
            if newid != record.id:
                print(f"{record.id} \n   changed to \n{newid} \n ---------- ",file=changed_records)

            record.id = newid
            record.description = newid
            SeqIO.write(record, outfile, "fasta")
            add_species_to_SH(output_species_dictionary, newid)



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
        

    modified_SH_data = []
    for sh in g_conflict_resolved:
        modified_SH_data.append({'SH': sh,
            'Old_Entries': ", ".join(input_species_dictionary[sh]),
            'New_Entries': ", ".join(output_species_dictionary[sh]),
            'Conflict_Type': "Different genera in SH"})

    for sh in s_conflict_resolved:
        modified_SH_data.append({'SH': sh,
            'Old_Entries': ", ".join(input_species_dictionary[sh]),
            'New_Entries': ", ".join(output_species_dictionary[sh]),
            'Conflict_Type': "Different species in SH"})

    for sh in s_and_g_conflict_resolved:
        modified_SH_data.append({'SH': sh,
            'Old_Entries': ", ".join(input_species_dictionary[sh]),
            'New_Entries': ", ".join(output_species_dictionary[sh]),
            'Conflict_Type': "Different species and genera in SH"})

    fieldnames = ['SH', 'Old_Entries', 'New_Entries', 'Conflict_Type']
    writer = csv.DictWriter(modified_SH, fieldnames = fieldnames)
    writer.writeheader()
    writer.writerows(modified_SH_data)



with open(summary_of_changes_txt, "w") as outfile:
    print(f"STATISTICS: Out of {statistics['total_records']} records, {statistics['all_levels_found_match']} have matches on GBIF on all taxonomical levels including species.", file=outfile)
    print(f"{statistics['s']} records had the species corrected by GBIF", file=outfile)
    print(f"{statistics['g']} records had the genus corrected by GBIF", file=outfile)
    print(f"{statistics['f']} records had the family corrected by GBIF", file=outfile)
    print(f"{statistics['o']} records had the order corrected by GBIF", file=outfile)
    print(f"{statistics['c']} records had the class corrected by GBIF", file=outfile)
    print(f"{statistics['p']} records had the phylum corrected by GBIF", file=outfile)
