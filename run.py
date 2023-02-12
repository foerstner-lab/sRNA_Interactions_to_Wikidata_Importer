from AnnotationToQIDTranslator import *
from RegulationClaimsImporter import *
print("*********** Welcome to the RegulonDB Interactions Wikidata importing tool ***********")
regulatory_file_path = input("Please enter the path of the RegulonDB interactions file: ")
GFF_file_path = input("Please enter the path of the GFF file: ")
Organism_QID = input("Please enter the QID of the organism: ")
RDB_obj = RegulonDB_Parser(regulatory_file_path)
parsed_interaction_df = RDB_obj.parse_file()
translated_obj = AnnotationToQIDTranslator(GFF_file_path, parsed_interaction_df, Organism_QID)
regulators_df, targets_df = translated_obj.build_translation_dictionary()
claims_df = translated_obj.prepare_claims(regulators_df, targets_df)
importer_obj = RegulationClaimsImporter(claims_df)
imports_count = importer_obj.import_claims()
print("Finished, " + str(imports_count) + " Claims imported")
