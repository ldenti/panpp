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

    seq_name = ""
    sfs_2 = set()
    seq_names_2 = set()

    for line in open(sfs_path_2):
        seq, s, l, *_ = line.strip("\n").split("\t")
        if seq != "*":
            seq_names_2.add(seq)
            if seq_name != "":
                if seq_name not in sfs_1:
                    print(f"{seq_name} not in SFS1")
                else:
                    if sfs_1[seq_name] != sfs_2:
                        print(seq_name, len(sfs_1[seq_name] & sfs_2))
                        print(seq_name, sfs_1[seq_name] - sfs_2)
                        print(seq_name, sfs_2 - sfs_1[seq_name])
                    else:
                        print(seq_name, "OK")
            seq_name = seq
            sfs_2 = set()
        sfs_2.add((s, l))

    seq_names_2.add(seq_name)
    if seq_name not in sfs_1:
        print(f"{seq_name} not in SFS1")
    else:
        if sfs_1[seq_name] != sfs_2:
            print(seq_name, len(sfs_1[seq_name] & sfs_2))
            print(seq_name, sfs_1[seq_name] - sfs_2)
            print(seq_name, sfs_2 - sfs_1[seq_name])
        else:
            print(seq_name, "OK")

    print(
        len(set(sfs_1.keys())) == len(seq_names_2)
        and len(seq_names_2) == len(set(sfs_1.keys()) & seq_names_2),
        len(set(sfs_1.keys())),
        len(seq_names_2),
        len(set(sfs_1.keys()) & seq_names_2),
    )
    for k in set(sfs_1.keys()) - seq_names_2:
        print(f"{k} in 1 but not in 2")
    for k in seq_names_2 - set(sfs_1.keys()):
        print(f"{k} in 2 but not in 1")


if __name__ == "__main__":
    main()
