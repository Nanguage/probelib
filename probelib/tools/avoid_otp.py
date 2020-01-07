import typing as t
from probelib.align.block import read_align_blocks, parse_region, count_overlap_with_region
from probelib.filters.avoid_otp import AvoidOTP


def avoid_otp(
        sam_path: str,
        out_path: str,
        target_regions: t.List[str],
        density_thresh: float = 1e-5,
        search_range: t.Tuple[int, int] = (-1e5, 1e5),
        ):
    regions = [parse_region(r) for r in target_regions]
    blocks = list(read_align_blocks(sam_path))

    # sort blocks
    def sort_key(b):
        return sum([count_overlap_with_region(b, r)[0] for r in regions])

    if regions:
        blocks.sort(key=sort_key, reverse=True)

    acc = AvoidOTP(regions, density_thresh, search_range)
    blocks = acc.filter(blocks)

    counted = []
    for b in blocks:
        if regions:
            c = [0, 0]
            for r in regions:
                c_ = count_overlap_with_region(b, r)
                c[0] += c_[0]
                c[1] += c_[1]
            if c[0] > 0:
                counted.append((b, c))
        else:
            c = (0, len(b[2]))
            counted.append((b, c))

    counted.sort(key=lambda t: t[1][0]/(t[1][0] + t[1][1]), reverse=True)
    with open(out_path, 'w') as f:
        for b, c in counted:
            name, seq, alns = b
            items = [name, seq, c[0], c[1], c[0]/sum(c)]
            items = [str(i) for i in items]
            outline = "\t".join(items) + "\n"
            f.write(outline)

    print(f"Done {out_path}")
