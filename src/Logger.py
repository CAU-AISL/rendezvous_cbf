import csv
import os
import numpy as np

class DataLogger:
    def __init__(self, log_dir="logs", file_name="log.csv"):
        self.log_dir = log_dir
        os.makedirs(log_dir, exist_ok=True)
        self.log_file_path = os.path.join(log_dir, file_name)
        self.header_written = False

    def _write_header(self, data_keys):
        with open(self.log_file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(data_keys)
        self.header_written = True

    def log(self, data):
        """
        Logs a dictionary of data to a CSV file.
        The keys of the dictionary will be used as the header if not already written.
        """
        data_keys = list(data.keys())
        if not self.header_written:
            self._write_header(data_keys)

        with open(self.log_file_path, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            row = [self._format_value(data[key]) for key in data_keys]
            writer.writerow(row)

    def _format_value(self, value):
        """
        Formats a value for CSV writing. Converts numpy arrays to space-separated strings.
        """
        if isinstance(value, np.ndarray):
            return ' '.join(map(str, value.flatten()))
        return str(value)

    def get_log_file_path(self):
        return self.log_file_path
