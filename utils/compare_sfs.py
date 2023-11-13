import sys


def parse_sfs(fpath):
    seq_name = ""
    specifics = {}
    for line in open(fpath):
        seq, s, l, *_ = line.strip("\n").split("\t")
        if seq != "*":
            seq_name = seq
            specifics[seq_name] = set()
        specifics[seq_name].add((s, l))
    return specifics


def main():
    sfs_path_1 = sys.argv[1]
    sfs_path_2 = sys.argv[2]

    sfs_1 = parse_sfs(sfs_path_1)
    sfs_2 = parse_sfs(sfs_path_2)

    assert sfs_1.keys() == sfs_2.keys()
    for k in sfs_1.keys():
        assert sfs_1[k] == sfs_2[k]
        # print(k, len(sfs_1[k] & sfs_2[k]))
        # print(sfs_1[k] - sfs_2[k])
        # print(sfs_2[k] - sfs_1[k])


if __name__ == "__main__":
    main()
