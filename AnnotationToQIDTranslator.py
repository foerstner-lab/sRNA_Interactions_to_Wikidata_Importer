from RegulonDB_Parser import *
from BCBio import GFF
from wikidataintegrator.wdi_core import WDItemEngine
from GFFRecordsMapper import *


class AnnotationToQIDTranslator:
    def __init__(self, GFF_file_path, parsed_interaction_df, Organism_QID):
        self.GFF_file_path = GFF_file_path
        self.parsed_interaction_df = parsed_interaction_df
        self.Organism_QID = Organism_QID

    def get_biotype(self, parent_rec, child_rec):
        parent_gbkey = ""
        child_gbkey = ""
        parent_biotype = ""
        child_ncrna_class = ""
        pseudo_flag = False
        bio_type = ""
        sub_types = []
        if 'gbkey' in parent_rec.qualifiers \
                and 'gbkey' in child_rec.qualifiers \
                and 'gene_biotype' in parent_rec.qualifiers:
            parent_gbkey = parent_rec.qualifiers['gbkey'][0]
            child_gbkey = child_rec.qualifiers['gbkey'][0]
            parent_biotype = parent_rec.qualifiers['gene_biotype'][0]
        else:
            print("Error getting the correct types of records, the file may not be supported")
            pass
        if 'ncrna_class' in child_rec.qualifiers:
            child_ncrna_class = child_rec.qualifiers['ncrna_class'][0]

        if parent_gbkey == "Gene":
            if 'pseudo' in parent_rec.qualifiers:
                if parent_rec.qualifiers['pseudo'][0] == 'true':
                    pseudo_flag = True
                else:
                    pseudo_flag = False
            elif 'pseudo' in parent_rec.qualifiers['gene_biotype'][0]:
                pseudo_flag = True
            else:
                pseudo_flag = False
            if 'protein_coding' == parent_biotype or 'CDS' == child_gbkey:
                bio_type = "protein"
                sub_types = ["protein"]
            elif 'RNA' in child_gbkey:
                if 'mRNA' == child_gbkey:
                    bio_type = "mRNA"
                    if parent_biotype == 'mRNA' or parent_biotype == 'mRNA_pseudogene':
                        sub_types = ['mRNA']
                    else:
                        print("Unknown sub type of mRNA")
                else:
                    bio_type = "ncRNA"
                    if 'rRNA' in child_gbkey:
                        if parent_biotype == 'rRNA' or parent_biotype == 'rRNA_pseudogene':
                            sub_types = ['rRNA']
                        else:
                            print("Unknown sub type of rRNA")
                    elif 'tRNA' in child_gbkey:
                        if parent_biotype == 'tRNA' or parent_biotype == 'tRNA_pseudogene':
                            sub_types = ['tRNA']
                        else:
                            print("Unknown sub type of tRNA")
                    elif 'tmRNA' in child_gbkey:
                        if parent_biotype == 'tmRNA' or parent_biotype == 'tmRNA_pseudogene':
                            sub_types = ['tmRNA']
                        else:
                            print("Unknown sub type of tmRNA")
                    elif 'ncRNA' in child_gbkey:
                        if parent_biotype == 'antisense_RNA' \
                                or child_ncrna_class == "antisense_RNA":
                            sub_types = ['antisense_RNA']
                        elif 'ncRNA' in parent_biotype \
                                and child_ncrna_class == 'other':
                            sub_types = ['unknown']
                        elif parent_biotype == 'SRP_RNA' \
                                and child_ncrna_class == 'SRP_RNA':
                            sub_types = ['SRP_RNA']
                        elif parent_biotype == 'RNase_P_RNA' \
                                and child_ncrna_class == 'RNase_P_RNA':
                            sub_types = ['RNase_P_RNA']
                    else:
                        print("Unsupported RNA type")
            else:
                print('error: cannot determine the bio type')
        else:
            print('error: cannot determine the bio type')
        if len(sub_types) == 0:
            print(child_rec.qualifiers)
        return bio_type, sub_types, pseudo_flag

    def build_query(self, params_dict):
        QUERY = ""
        query_file = open('query_templates/FIND_QID_QUERY_HEADER.rq', 'r')
        query_template = query_file.read()
        query_file.close()
        QUERY = query_template
        QUERY = QUERY.replace("#QID#", self.Organism_QID)
        QUERY = QUERY.replace("#LOCUS_TAG#", params_dict['locus_tag'])
        QUERY = QUERY.replace("#GENE_NAME#", params_dict['gene_name'])
        QUERY += "\n{"
        if 'RNA' in params_dict['biotype']:
            QUERY += "\n{?item wdt:P31 wd:Q11053.}" # instance of / subclass of RNA
            QUERY += "\nUNION\n{?item wdt:P279 wd:Q11053.}"  # instance of / subclass of RNA
            if params_dict['biotype'] == 'ncRNA':
                QUERY += "\nUNION\n{?item wdt:P31 wd:Q427087.}"  # instance / subclass of of non-coding RNA
                QUERY += "\nUNION\n{?item wdt:P279 wd:Q427087.}"  # instance / subclass of of non-coding RNA
                for type in params_dict['sub_types']:
                    if type == 'tRNA':
                        QUERY += "\nUNION\n{?item wdt:P31 wd:Q201448.}"  # instance of / subclass of tRNA
                        QUERY += "\nUNION\n{?item wdt:P279 wd:Q201448.}"  # instance of / subclass of tRNA
                    elif type == 'tmRNA':
                        QUERY += "\nUNION\n{?item wdt:P31 wd:Q285904.}"  # instance of / subclass of tmRNA
                        QUERY += "\nUNION\n{?item wdt:P279 wd:Q285904.}"  # instance of / subclass of tmRNA
                    elif type == 'antisense_RNA':
                        QUERY += "\nUNION\n{?item wdt:P31 wd:Q423832.}"  # instance of / subclass of antisense RNA
                        QUERY += "\nUNION\n{?item wdt:P279 wd:Q423832.}"  # instance of / subclass of antisense RNA
                    elif type == 'unknown':
                        QUERY += "\nUNION\n{?item wdt:P31 wd:Q24238356.}"  # instance of / subclass of unknown RNA
                        QUERY += "\nUNION\n{?item wdt:P279 wd:Q24238356.}"  # instance of / subclass of unknown RNA
                    elif type == 'SRP_RNA':
                        QUERY += "\nUNION\n{?item wdt:P31 wd:Q424665.}"  # instance of / subclass of signal recognition particle
                        QUERY += "\nUNION\n{?item wdt:P279 wd:Q424665.}"  # instance of / subclass of signal recognition particle
                    elif type == 'RNase_P_RNA':
                        QUERY += "\nUNION\n{?item wdt:P31 wd:Q1012651.}"  # instance of / subclass of RNase P
                        QUERY += "\nUNION\n{?item wdt:P279 wd:Q1012651.}"  # instance of / subclass of RNase P
                    else:
                        print("Unsupported type of RNA")
            elif params_dict['biotype'] == 'mRNA':
                QUERY += "\n{?item wdt:P31 wd:Q188928.}"  # instance of / subclass of mRNA
                QUERY += "\nUNION\n{?item wdt:P279 wd:Q188928.}"  # instance of / subclass of mRNA
            else:
                print("Un-recognized subtype of RNA.")
        elif params_dict['biotype'] == 'gene':
            QUERY += "\n{?item wdt:P31 wd:Q7187.}"  # instance of / subclass of gene
            QUERY += "\nUNION\n{?item wdt:P279 wd:Q7187.}"  # instance of / subclass of gene
            QUERY += "\nUNION\n{?item wdt:P31 wd:Q20747295.}"  # instance of / subclass of protein coding gene
            QUERY += "\nUNION\n{?item wdt:P279 wd:Q20747295.}"  # instance of / subclass of protein coding gene
            QUERY += "\nUNION\n{?item wdt:P31 wd:Q277338.}"  # instance of / subclass of pseudogene
            QUERY += "\nUNION\n{?item wdt:P279 wd:Q277338.}"  # instance of / subclass of pseudogene
        elif params_dict['biotype'] == 'protein':
            QUERY += "\n{?item wdt:P31 wd:Q8054.}"  # instance of / subclass of protein
            QUERY += "\n{?item wdt:P279 wd:Q8054.}"  # instance of / subclass of protein
        else:
            print("Un-recognized type.")
        QUERY += "\n}"

        QUERY += '\n{\n{?item wdt:P2393 "' + params_dict['locus_tag'] + '".}'  # NCBI Locus Tag
        if 'GeneID' in params_dict['Dbxref']:
            QUERY += '\nUNION\n{?item wdt:P351 "' + params_dict['Dbxref']['GeneID'] + '".}'  # Entrez Gene ID
        elif 'Genbank' in params_dict['Dbxref']:
            QUERY += '\nUNION\n{?item wdt:P637 "' + params_dict['Dbxref']['Genbank'] + '".}'  # RefSeq Protein ID
        elif 'UniProtKB/Swiss-Prot' in params_dict['Dbxref']:
            QUERY += '\nUNION\n{?item wdt:P352 "' + params_dict['Dbxref']['UniProtKB/Swiss-Prot'] + '".}'  # UniProt protein ID
        # elif 'ASAP' in params_dict['Dbxref']:
        #     QUERY += '\nUNION\n{?item wdt:P?? "' + params_dict['Dbxref']['ASAP'] + '".}'  # ??
        # elif 'EcoGene' in params_dict['Dbxref']:
        #     QUERY += '\nUNION\n{?item wdt:P?? "' + params_dict['Dbxref']['EcoGene'] + '".}'  # ??
        else:
            print("Unsupported type of Identifiers")
        QUERY += "\n}\n}"
        return QUERY

    def get_QID(self, query):
        item_QID = []
        results = WDItemEngine.execute_sparql_query(query)['results']['bindings']
        if len(results) == 0:
            item_QID.append("NOT_FOUND_IN_WD")
        elif len(results) == 1:
            for result in results:
                item_QID.append(result['item']['value'].replace("http://www.wikidata.org/entity/", ""))
        else:
            for result in results:
                item_QID.append(result['item']['value'].replace("http://www.wikidata.org/entity/", ""))
                print("Warning: Query returns more than one item for the same gene name! Selected: "
                      + result['item']['value'])
        return item_QID

    def build_translation_dictionary(self):

        regulator_list = []
        target_list = []
        parsed_gff_obj = GFF.parse(open(self.GFF_file_path, 'r'), target_lines=1)
        for index, row in self.parsed_interaction_df.iterrows():
            if row['regulator'] not in regulator_list:
                regulator_list.append(row['regulator'])
            if row['target'] not in target_list:
                target_list.append(row['target'])
        regulators_df = pd.DataFrame({'gene_name': regulator_list, 'QID': None})
        targets_df = pd.DataFrame({'gene_name': target_list, 'QID': None})
        del regulator_list
        del target_list
        GFFRecordsMapper_obj = GFFRecordsMapper(parsed_gff_obj)
        records_df = GFFRecordsMapper_obj.map_gene_data()
        for regulator_index, regulator in regulators_df.iterrows():
            query = ""
            GFF_found = False
            Dbxref_tmp = {}
            for index, row in records_df.iterrows():
                if regulator['gene_name'] == row['Parent_Record'].qualifiers['Name'][0] \
                        or regulator['gene_name'] == row['Parent_Record'].qualifiers['gene'][0] \
                        or regulator['gene_name'] in row['Parent_Record'].qualifiers['gene_synonym']:
                    GFF_found = True
                    for db_id in row['Record'].qualifiers['Dbxref']:
                        if 'GeneID' in db_id:
                            Dbxref_tmp['GeneID'] = db_id.split(':', 1)[-1]
                        elif 'Genbank' in db_id:
                            Dbxref_tmp['Genbank'] = db_id.split(':', 1)[-1]
                        elif 'UniProtKB/Swiss-Prot' in db_id:
                            Dbxref_tmp['UniProtKB/Swiss-Prot'] = db_id.split(':', 1)[-1]
                        elif 'ASAP' in db_id:
                            Dbxref_tmp['ASAP'] = db_id.split(':', 1)[-1]
                        elif 'EcoGene' in db_id:
                            Dbxref_tmp['ASAP'] = db_id.split(':', 1)[-1]
                        else:
                            print("Unsupported type of identifiers " + db_id)
                    biotype, sub_types, pseudo = self.get_biotype(row['Parent_Record'], row['Record'])
                    query = self.build_query({'gene_name': regulator['gene_name'], 'biotype': biotype,
                                              'sub_types': sub_types,
                                              'locus_tag': row['Parent_Record'].qualifiers['locus_tag'][0],
                                              'Dbxref': Dbxref_tmp})
                    break
            if not GFF_found:
                regulators_df.iloc[regulator_index]['QID'] = 'NOT_FOUND_IN_GFF'
            else:
                regulators_df.iloc[regulator_index]['QID'] = self.get_QID(query)
        for target_index, target in targets_df.iterrows():
            query = ""
            GFF_found = False
            Dbxref_tmp = {}
            for index, row in records_df.iterrows():
                if target['gene_name'] == row['Parent_Record'].qualifiers['Name'][0] \
                        or target['gene_name'] == row['Parent_Record'].qualifiers['gene'][0] \
                        or target['gene_name'] in row['Parent_Record'].qualifiers['gene_synonym']:
                    GFF_found = True
                    for db_id in row['Parent_Record'].qualifiers['Dbxref']:
                        if 'GeneID' in db_id:
                            Dbxref_tmp['GeneID'] = db_id.split(':', 1)[-1]
                        elif 'Genbank' in db_id:
                            Dbxref_tmp['Genbank'] = db_id.split(':', 1)[-1]
                        elif 'UniProtKB/Swiss-Prot' in db_id:
                            Dbxref_tmp['UniProtKB/Swiss-Prot'] = db_id.split(':', 1)[-1]
                        elif 'ASAP' in db_id:
                            Dbxref_tmp['ASAP'] = db_id.split(':', 1)[-1]
                        elif 'EcoGene' in db_id:
                            Dbxref_tmp['ASAP'] = db_id.split(':', 1)[-1]
                        else:
                            print("Unsupported type of identifiers " + db_id)
                    # Target must be a gene record
                    # biotype, sub_types, pseudo = self.get_biotype(row['Parent_Record'], row['Record'])
                    query = self.build_query({'gene_name': target['gene_name'], 'biotype': 'gene',
                                              'sub_types': ['gene'],
                                              'locus_tag': row['Parent_Record'].qualifiers['locus_tag'][0],
                                              'Dbxref': Dbxref_tmp})
                    break
            if not GFF_found:
                targets_df.iloc[target_index]['QID'] = 'NOT_FOUND_IN_GFF'
            else:
                targets_df.iloc[target_index]['QID'] = self.get_QID(query)
        return regulators_df, targets_df

    def set_property(self, regulation_type):
        r_property = ""
        if regulation_type == "antisense":
            r_property = "P3777"  # antisense inhibitor
        elif regulation_type == "" or regulation_type == 'unknown' or regulation_type is None:
            r_property = "P128"  # Regulates
        elif regulation_type == 'represses processing':
            r_property = "P128"  # Regulates
        elif regulation_type == 'inactivated':
            r_property = "P128"  # Regulates
        elif regulation_type == 'base-pairing':
            r_property = "P128"  # Regulates
        elif regulation_type == 'translocation':
            r_property = "P128"  # Regulates
        elif regulation_type == 'catalytic part of RNase P':
            r_property = "P128"  # Regulates
        elif regulation_type == 'stability of the transcript':
            r_property = "P128"  # Regulates
        elif regulation_type == "antagonist":
            r_property = "P3773"  # antagonist
        elif regulation_type == "activated":
            r_property = "P3771"  # activator of
        else:
            r_property = ""
            print("Warning: Un-recognized type of regulation of: '" + regulation_type + "'")
        return r_property

    def prepare_claims(self, regulators_df, targets_df):
        tmp_list_dict = []
        for index, row in self.parsed_interaction_df.iterrows():
            tmp_dict = {'regulator_QID': "", 'property': "", 'target_QID': ""}
            for reg_index, regulator in regulators_df.iterrows():
                if regulator['QID'] != "NOT_FOUND_IN_GFF":
                    if regulator['gene_name'] == row['regulator']:
                        if regulator['QID'][0] != "NOT_FOUND_IN_WD":
                            tmp_dict['regulator_QID'] = regulator['QID'][0]
                        else:
                            tmp_dict['regulator_QID'] = ""
                        break
            for tar_index, target in targets_df.iterrows():
                if target['QID'] != "NOT_FOUND_IN_GFF":
                    if target['gene_name'] == row['target']:
                        if target['QID'][0] != "NOT_FOUND_IN_WD":
                            for QID in target['QID']:
                                tmp_dict['target_QID'] = QID
                                tmp_dict['property'] = self.set_property(row['regulation-type'])
                                if tmp_dict['regulator_QID'] != ""\
                                        and tmp_dict['property'] != ""\
                                        and tmp_dict['target_QID'] != "":
                                    tmp_list_dict.append(tmp_dict)
                        break
            if tmp_dict['regulator_QID'] == "":
                print("Warning: Gene '" + row['regulator'] + "' is not found in the GFF file or Wikidata and ignored.")
            if tmp_dict['target_QID'] == "":
                print("Warning: Gene '" + row['target'] + "' is not found in the GFF file or Wikidata and ignored.")
        claims_df = pd.DataFrame.from_records(tmp_list_dict)
        return claims_df

