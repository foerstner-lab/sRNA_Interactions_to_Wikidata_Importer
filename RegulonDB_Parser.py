import pandas as pd
import io
import csv

class RegulonDB_Parser:

    def __init__(self, regulatory_file_path):
        self.regulatory_file_path = regulatory_file_path
        self.total_lines_count = 0
        self.ignored_lines_count = 0
        self.filtered_lines_count = 0

    def parse_file(self):

        file_filtered_content = "regulonDB-id\tregulator\ttarget\tsite-left\tsite-right\ttarget-strand\tfunction\t" \
                                "bind-site-seq\tregulation-type\tmechanism\tevidence\tconfidence-level\n"
        regulatory_file_obj = open(self.regulatory_file_path, 'r')
        for line in regulatory_file_obj:
            self.total_lines_count += 1
            if line.startswith('#'):
                self.ignored_lines_count += 1
            else:
                file_filtered_content += line
                self.filtered_lines_count += 1
        return pd.DataFrame.from_records(list(csv.DictReader(io.StringIO(file_filtered_content), delimiter='\t')))

    def displayDataCounted(self):
        print("Total lines: " + str(self.total_lines_count) + ", header lines: " + str(
            self.ignored_lines_count) + ", data lines:" + str(self.filtered_lines_count))