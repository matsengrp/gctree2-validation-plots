import os
import python.paircluster as paircluster
import argparse

parser = argparse.ArgumentParser(
    description="Produce gctree inputs from a Partis output directory"
)
parser.add_argument(
    "input_dir", type=str, help="The partis output directory containing simulations"
)
parser.add_argument(
    "output_dir",
    type=str,
    help="The directory in which to place gctree inputs, in one subdirectory per simulation",
)
args = parser.parse_args()

pc = paircluster.read_paired_dir(args.input_dir)
annotations_h = pc[("igh", "igk")]["antn_lists"]["igh"]
annotations_k = pc[("igh", "igk")]["antn_lists"]["igk"]

for simu_idx, (annotation_h, annotation_k) in enumerate(
    zip(annotations_h, annotations_k)
):
    this_simu_dir = args.output_dir + "/" + str(simu_idx) + "/"
    os.makedirs(this_simu_dir)
    naive_h = annotation_h["naive_seq"]
    igh_length = len(naive_h)
    naive_concat = naive_h + annotation_k["naive_seq"]

    frame_igh = (annotation_h["v_5p_del"] + 1) % 3
    frame_igk = (annotation_k["v_5p_del"] + 1) % 3

    k_seq_alignment = {
        uid: seq for uid, seq in zip(annotation_k["unique_ids"], annotation_k["seqs"])
    }

    alignment_concat = {}
    for seqid_h, paired_k_seqids, seq_h in zip(
        annotation_h["unique_ids"], annotation_h["paired-uids"], annotation_h["seqs"]
    ):
        alignment_concat[seqid_h[:-4]] = seq_h + k_seq_alignment[paired_k_seqids[0]]

    with open(this_simu_dir + "frame1_igh.txt", "w") as fh:
        fh.write(str(frame_igh))

    with open(this_simu_dir + "frame2_igk.txt", "w") as fh:
        fh.write(str(frame_igk))

    with open(this_simu_dir + "chain_split_igh_len.txt", "w") as fh:
        fh.write(str(len(naive_h)))

    with open(this_simu_dir + "all_seqs.fasta", "w") as fh:
        for seqid, seq in alignment_concat.items():
            fh.write(">" + seqid + "\n" + seq + "\n")

    with open(this_simu_dir + "inference.fasta", "w") as fh:
        fh.write(">naive\n" + naive_concat + "\n")
        for seqid, seq in alignment_concat.items():
            if "leaf" in seqid:
                fh.write(">" + seqid + "\n" + seq + "\n")

    with open(this_simu_dir + "true_tree.nwk", "w") as fh:
        fh.write(annotation_h["tree"].replace("-igh", "").replace("[&R] ", ""))
