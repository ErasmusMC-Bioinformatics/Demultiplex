import argparse
from collections import defaultdict
import logging
import sys
import os
from pprint import pprint
from fuzzywuzzy import process
from Bio.Seq import Seq
from Bio import SeqIO

MIDS = {
    "MID1": "ACGAGTGCGT",
    "MID2": "ACGCTCGACA",
    "MID3": "AGACGCACTC",
    "MID4": "AGCACTGTAG",
    "MID5": "ATCAGACACG",
    "MID6": "ATATCGCGAG",
    "MID7": "CGTGTCTCTA",
    "MID8": "CTCGCGTGTC",
    "MID10": "TCTCTATGCG",
    "MID11": "TGATACGTCT",
    "MID13": "CATAGTAGTG",
    "MID14": "CGAGAGATAC",
    "MID15": "ATACGACGTA",
    "MID16": "TCACGTACTA",
    "MID17": "CGTCTAGTAC",
    "MID18": "TCTACGTAGC",
    "MID19": "TGTACTACTC",
    "MID20": "ACGACTACAG",
    "MID21": "CGTAGACTAG",
    "MID22": "TACGAGTATG",
    "MID23": "TACTCTCGTG",
    "MID24": "TAGAGACGAG",
    "MID25": "TCGTCGCTCG",
    "MID26": "ACATACGCGT",
    "MID27": "ACGCGAGTAT",
    "MID28": "ACTACTATGT",
    "MID29": "ACTGTACAGT",
    "MID30": "AGACTATACT",
    "MID31": "AGCGTCGTCT",
    "MID32": "AGTACGCTAT",
    "MID33": "ATAGAGTACT",
    "MID34": "CACGCTACGT",
    "MID35": "CAGTAGACGT",
    "MID36": "CGACGTGACT",
    "MID37": "TACACACACT",
    "MID38": "TACACGTGAT",
    "MID39": "TACAGATCGT",
    "MID40": "TACGCTGTCT",
    "MID41": "TAGTGTAGAT",
    "MID42": "TCGATCACGT",
    "MID43": "TCGCACTAGT",
    "MID44": "TCTAGCGACT",
    "MID45": "TCTATACTAT",
    "MID46": "TGACGTATGT",
    "MID47": "TGTGAGTAGT",
    "MID48": "ACAGTATATA",
    "MID49": "ACGCGATCGA",
    "MID50": "ACTAGCAGTA",
    "MID51": "AGCTCACGTA",
    "MID52": "AGTATACATA",
    "MID53": "AGTCGAGAGA",
    "MID54": "AGTGCTACGA",
    "MID55": "CGATCGTATA",
    "MID56": "CGCAGTACGA",
    "MID57": "CGCGTATACA",
    "MID58": "CGTACAGTCA",
    "MID59": "CGTACTCAGA",
    "MID60": "CTACGCTCTA",
    "MID61": "CTATAGCGTA",
    "MID62": "TACGTCATCA",
    "MID63": "TAGTCGCATA",
    "MID64": "TATATATACA",
    "MID65": "TATGCTAGTA",
    "MID66": "TCACGCGAGA",
    "MID67": "TCGATAGTGA",
    "MID68": "TCGCTGCGTA",
    "MID69": "TCTGACGTCA",
    "MID70": "TGAGTCAGTA",
    "MID71": "TGTAGTGTGA",
    "MID72": "TGTCACACGA",
    "MID73": "TGTCGTCGCA",
    "MID74": "ACACATACGC",
    "MID75": "ACAGTCGTGC",
    "MID76": "ACATGACGAC",
    "MID77": "ACGACAGCTC",
    "MID78": "ACGTCTCATC",
    "MID79": "ACTCATCTAC"
}

def read_fastq(fastq_path):
    with open(fastq_path, 'rU') as handle:
        name, seq, thing, qual = None, None, "+", None
        for record in SeqIO.parse(handle, "fastq"):
            name = record.id
            seq = record.seq
            qual = record.letter_annotations
            yield (name, seq, thing, qual)

def main():
    logging.basicConfig(filename="./log.html", level=logging.INFO, format="%(asctime)s:&emsp;%(message)s <br />", datefmt='%Y/%m/%d %H:%M:%S')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    logging.info("Started demultiplexing")

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fastq", help="")
    parser.add_argument("--output-dir", help="")
    parser.add_argument("--mid-file", help="")

    args = parser.parse_args()

    input_fastq = args.input_fastq
    bn_fastq = os.path.basename(input_fastq)
    output_dir = args.output_dir
    mid_file = args.mid_file

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fastq_dir = os.path.join(output_dir, "fastq")
    fasta_dir = os.path.join(output_dir, "fasta")

    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)

    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)

    logging.info("Input: {0}".format(input_fastq))

    output_dict = {}

    with open(mid_file, 'r') as m: # create a dict that, for every sample, holds the mids and output files.
        for mid_line in m:
            splt = mid_line.split("\t")
            sample = splt[0]
            mid = splt[1].rstrip()

            output_fastq = os.path.join(fastq_dir, "{0}.fastq".format(sample))
            output_fasta = os.path.join(fasta_dir, "{0}.fasta".format(sample))
            output_disc = os.path.join(fastq_dir, "{0}_discarded.fastq".format(os.path.splitext(bn_fastq)[0]))

            output_dict[sample] = {
                "tag": MIDS[mid],
                "rev": str(Seq(MIDS[mid]).reverse_complement()),
                "mid": mid,
                "output_fastq": output_fastq,
                "output_fasta": output_fastq,
                "fastq": open(output_fastq, 'w'),
                "fasta": open(output_fasta, 'w'),
                "discarded": open(output_disc, 'w')
            }

    total_seq_count = 0
    total_assigned = 0
    per_sample_count_dict = defaultdict(int)

    fuzzy_forward_matches = 0
    IDs = []
    found_location = defaultdict(int)
    for ID, seq, thing, qual in read_fastq(input_fastq):
        ID = ID[1:]
        seq_half = int(len(seq) / 2)
        total_seq_count += 1

        for sample, dic in output_dict.items():
            if seq.find(dic["tag"], 0, seq_half) != -1:
                per_sample_count_dict[sample] += 1
                found_location[str(seq.find(dic["tag"]))] += 1
                per_sample_count_dict["{0}_forward_fh".format(sample)] += 1
                dic["fastq"].write("@{0}\n{1}\n{2}\n{3}\n".format(ID, seq, thing, qual))
                dic["fasta"].write(">{0}\n{1}\n".format(ID, seq))
                total_assigned += 1
                IDs.append("{0} [('{1}', '{2}', 90, 'f_f45')] exact".format(ID, sample, dic['tag']))
                break
            if seq.find(dic["rev"], seq_half, len(seq)) != -1:
                per_sample_count_dict[sample] += 1
                per_sample_count_dict["{0}_reverse_sh".format(sample)] += 1
                dic["fastq"].write("@{0}\n{1}\n{2}\n{3}\n".format(ID, seq, thing, qual))
                dic["fasta"].write(">{0}\n{1}\n".format(ID, seq))
                total_assigned += 1
                IDs.append("{0} [('{1}', '{2}', 90, 'r_r45')] exact".format(ID, sample, dic['rev']))
                break
            if seq.find(dic["rev"], 0, seq_half) != -1:
                per_sample_count_dict[sample] += 1
                per_sample_count_dict["{0}_reverse_fh".format(sample)] += 1
                dic["fastq"].write("@{0}\n{1}\n{2}\n{3}\n".format(ID, seq, thing, qual))
                dic["fasta"].write(">{0}\n{1}\n".format(ID, seq))
                total_assigned += 1
                IDs.append("{0} [('{1}', '{2}', 90, 'r_f45')] exact".format(ID, sample, dic['rev']))
                break
            if seq.find(dic["tag"], seq_half, len(seq)) != -1:
                per_sample_count_dict[sample] += 1
                per_sample_count_dict["{0}_forward_sh".format(sample)] += 1
                dic["fastq"].write("@{0}\n{1}\n{2}\n{3}\n".format(ID, seq, thing, qual))
                dic["fasta"].write(">{0}\n{1}\n".format(ID, seq))
                total_assigned += 1
                IDs.append("{0} [('{1}', '{2}', 90, 'f_r45')] exact".format(ID, sample, dic['tag']))
                break

        else:  # only if there wasn't an exact match, no 'break'
            dic["discarded"].write("@{0}\n{1}\n{2}\n{3}\n".format(ID, seq, thing, qual))

        if (total_seq_count % 100000) == 0:
            logging.info("Did {0}, {1} assigned with {2} fuzzy assigned".format(total_seq_count, total_assigned, fuzzy_forward_matches))

    log_file = "{0}_log.txt".format(os.path.splitext(bn_fastq)[0])
    log_file = os.path.join(output_dir, log_file)
    total_f_mid_fh = 0
    total_r_mid_sh = 0
    total_f_mid_sh = 0
    total_r_mid_fh = 0
    with open(log_file, 'w') as lf:
        lf.write("{0}\n".format("\t".join(["input", "sample", "assigned", "percentage", "f_mid_fh", "r_mid_sh", "f_mid_sh", "r_mid_fh"])))
        for sample, dic in output_dict.items():
            dic["fastq"].close()
            dic["fasta"].close()
            f_mid_fh = per_sample_count_dict["{0}_forward_fh".format(sample)]
            r_mid_sh = per_sample_count_dict["{0}_reverse_sh".format(sample)]
            f_mid_sh = per_sample_count_dict["{0}_forward_sh".format(sample)]
            r_mid_fh = per_sample_count_dict["{0}_reverse_fh".format(sample)]

            total_f_mid_fh += f_mid_fh
            total_r_mid_sh += r_mid_sh
            total_f_mid_sh += f_mid_sh
            total_r_mid_fh += r_mid_fh
            lf.write("{0}\n".format("\t".join([str(total_seq_count), sample, str(per_sample_count_dict[sample]), str(round(per_sample_count_dict[sample] / float(total_seq_count) * 100.0, 3)), str(str(f_mid_fh)), str(str(r_mid_sh)), str(str(f_mid_sh)), str(str(r_mid_fh))])))
        lf.write("{0}\n".format("\t".join([str(total_seq_count), "total", str(total_assigned), str(round(total_assigned / float(total_seq_count) * 100.0, 3)), str(str(total_f_mid_fh)), str(str(total_r_mid_sh)), str(str(total_f_mid_sh)), str(str(total_r_mid_fh))])))

    print(total_seq_count)
    print(total_assigned)
    print("{0}%".format(round(total_assigned / float(total_seq_count) * 100.0)))
    pprint(dict(per_sample_count_dict))
    #pprint(dict(found_location))
    """
    with open("pauline_forward_exact_ids.txt", 'w') as o:
        for ID in IDs:
            o.write("{0}\n".format(ID))
    """

if __name__ == "__main__":
    main()
