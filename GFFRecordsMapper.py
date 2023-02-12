import pandas as pd
from itertools import cycle


class GFFRecordsMapper:

    def __init__(self, gff_parsed_obj):
        self.gff_parsed_obj = gff_parsed_obj

    def map_gene_data(self):
        cnt = 0
        gff_lst = list(self.gff_parsed_obj)
        gff_len = len(gff_lst)
        records_mapped = []
        for record in cycle(gff_lst):
            cnt += 1
            if gff_len*2+1 == cnt:
                break
            if cnt <= gff_len:
                for feature in record.features:

                    if feature.qualifiers['gbkey'][0] != "Src" \
                            and feature.qualifiers['gbkey'][0] != "repeat_region" \
                            and feature.qualifiers['gbkey'][0] != "mobile_element" \
                            and feature.qualifiers['gbkey'][0] != "misc_feature" \
                            and feature.qualifiers['gbkey'][0] != "rep_origin" \
                            and feature.type != "exon":
                        if 'Parent' in feature.qualifiers:
                            records_mapped.append({'GFF_ID': record.id, 'Record': feature})
            if cnt > gff_len: # Second cycle
                for item in records_mapped:
                    if 'Parent_Record' not in item:
                        for feat in record.features:
                            if feat.qualifiers['gbkey'][0] != "Src" \
                                    and feat.qualifiers['gbkey'][0] != "repeat_region" \
                                    and feat.qualifiers['gbkey'][0] != "mobile_element" \
                                    and feat.qualifiers['gbkey'][0] != "misc_feature" \
                                    and feat.qualifiers['gbkey'][0] != "rep_origin" \
                                    and feat.type != "exon"\
                                    and 'Parent' not in feat.qualifiers:
                                if feat.qualifiers['ID'][0] == item['Record'].qualifiers['Parent'][0]:
                                    if item['Record'].location == feat.location:
                                        item['Parent_Record'] = feat
                                    elif item['Record'].location.start == feat.location.start \
                                            and item['Record'].location.end != feat.location.end:
                                        item['Parent_Record'] = feat
                                    elif item['Record'].location.end == feat.location.end \
                                            and item['Record'].location.start != feat.location.start:
                                        item['Parent_Record'] = feat
                                    elif item['Record'].location.start > feat.location.start \
                                            and item['Record'].location.end < feat.location.end:
                                        item['Parent_Record'] = feat
        Mapped_df = pd.DataFrame.from_dict(records_mapped)
        return Mapped_df
