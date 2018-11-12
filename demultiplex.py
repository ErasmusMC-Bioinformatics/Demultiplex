import argparse
import csv
import logging
import os
import sys
from collections import defaultdict
from collections import namedtuple

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def sniff_format(file_path):
    """
    Try to guess the file format (fastq/fasta) by looking at the first character of the first line.
    Should be '@' for fastq and '>' for fasta.
    """
    with open(file_path, 'r') as file_handle:
        for line in file_handle:
            if line.startswith("@"):
                return "fastq"
            if line.startswith(">"):
                return "fasta"
            break
        return None


def search_barcode_in_first_half(sequence, barcode):
    if type(sequence) is Seq:
        sequence = str(sequence).lower()
    elif type(sequence) is SeqRecord:
        sequence = str(sequence.seq).lower()
    return sequence.find(barcode, 0, int(len(sequence) / 2))


def search_barcode_in_second_half(sequence, barcode):
    if type(sequence) is Seq:
        sequence = str(sequence).lower()
    elif type(sequence) is SeqRecord:
        sequence = str(sequence.seq).lower()
    return sequence.find(barcode, int(len(sequence) / 2))


def search_barcodes_in_sequence(barcode_datas, sequence):
    for barcode_data in barcode_datas:
        barcode_search_position = search_barcode_in_first_half(sequence, barcode_data.barcode)
        if barcode_search_position != -1:
            return barcode_search_position, barcode_data, False
        barcode_search_position = search_barcode_in_second_half(sequence, barcode_data.barcode_reverse)
        if barcode_search_position != -1:
            return barcode_search_position, barcode_data, True

        barcode_search_position = search_barcode_in_first_half(sequence, barcode_data.barcode_reverse)
        if barcode_search_position != -1:
            return barcode_search_position, barcode_data, True

        barcode_search_position = search_barcode_in_second_half(sequence, barcode_data.barcode)
        if barcode_search_position != -1:
            return barcode_search_position, barcode_data, False

    return -1, None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="The input file", required=True)
    parser.add_argument("-f", "--format", help="The format of the input file (fastq/fasta)", default="auto",
                        choices=["fasta", "fastq", "auto"])
    parser.add_argument("-o", "--output-dir", help="The output dir", required=True)
    parser.add_argument("-m", "--mapping-file",
                        help="A tab seperated file containing two columns, ID and barcode (no header)", required=True)

    args = parser.parse_args()

    input_file_path = args.input
    basename_input_file_path = os.path.basename(input_file_path)
    input_basename_no_ext, input_extension = os.path.splitext(basename_input_file_path)

    input_format = args.format

    logging.basicConfig(stream=sys.stdout, level=logging.INFO)

    if input_format == "auto":
        if input_extension in [".fasta", ".fa"]:
            input_format = "fasta"
        elif input_extension in [".fastq", ".fq"]:
            input_format = "fastq"
        else:
            logging.info("Can't auto detect input format based on extension: {0}".format(input_extension))
            logging.info("Sniffing format...")
            input_format = sniff_format(input_file_path)
            if not input_format:
                logging.error("Can't auto detect input format")
                sys.exit(1)
            logging.info("Sniffed '{0}' as format.".format(input_format))

    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(args.mapping_file, newline="") as handle:
        ID_barcode_mapping = list(csv.DictReader(handle, fieldnames=['ID', 'barcode'], delimiter="\t"))

    logging.info("Input:   {0}".format(input_file_path))
    logging.info("Format:  {0}".format(input_format))
    logging.info("Output:  {0}".format(output_dir))
    logging.info("Mapping: {0} ({1} mappings)".format(args.mapping_file, len(ID_barcode_mapping)))

    BarcodeData = namedtuple("BarcodeData", [
        "ID",
        "barcode",
        "barcode_reverse",
        "output_file_path",
        "output_file_handle"
    ])

    barcode_data_dict = defaultdict(list)
    ID_file_handle_dict = {}

    for ID_barcode in ID_barcode_mapping:
        ID = ID_barcode["ID"]
        barcode = ID_barcode["barcode"]

        logging.info("{0}:\t\t{1}".format(ID, barcode))

        output_file_path = os.path.join(
            output_dir,
            "{0}.{1}".format(ID, input_format)
        )

        if ID not in ID_file_handle_dict:
            ID_file_handle = open(output_file_path, 'w')
            ID_file_handle_dict[ID] = ID_file_handle

        ID_file_handle = ID_file_handle_dict[ID]

        barcode_data_dict[ID] += [BarcodeData(
            ID=ID,
            barcode=barcode.lower(),
            barcode_reverse=str(Seq(barcode, generic_dna).reverse_complement()).lower(),
            output_file_path=output_file_path,
            output_file_handle=ID_file_handle
        )]

    discarded_output_file_path = os.path.join(
        output_dir,
        "{0}.{1}".format("discarded", input_format)
    )

    total_sequences = 0
    sequences_assigned_by_id = defaultdict(int)

    with open(input_file_path, 'r') as input_file_handle, open(discarded_output_file_path,
                                                               'w') as discarded_output_handle:
        for record in SeqIO.parse(input_file_handle, input_format):
            total_sequences += 1
            for ID, barcode_datas in barcode_data_dict.items():
                barcode_position, barcode_data, reverse = search_barcodes_in_sequence(barcode_datas, record)
                if barcode_position == -1:
                    continue
                barcode = barcode_data.barcode if not reverse else barcode_data.barcode_reverse
                sequences_assigned_by_id[barcode_data.ID] += 1
                logging.debug(str(record.seq))
                logging.debug(" " * barcode_position + barcode)
                SeqIO.write(record, ID_file_handle_dict[ID], input_format)
                break
            else:  # no match # TODO fuzzy match ?
                SeqIO.write(record, discarded_output_handle, input_format)
            if total_sequences % 10000 == 0:
                assigned_sequences = sum(sequences_assigned_by_id.values())
                logging.info("Processed {0} sequences, assigned {1} ({2}%)".format(
                    total_sequences,
                    assigned_sequences,
                    round(assigned_sequences / total_sequences * 100)
                ))

    for file_handle in ID_file_handle_dict.values():
        file_handle.close()

    assigned_sequences = sum(sequences_assigned_by_id.values())
    logging.info("Processed {0} sequences, assigned {1} ({2}%)".format(
        total_sequences,
        assigned_sequences,
        round(assigned_sequences / total_sequences * 100)
    ))


if __name__ == "__main__":
    main()

