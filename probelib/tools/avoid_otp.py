import typing as t
from probelib.align.block import read_align_blocks, parse_region, count_ovlap_with_region
from probelib.filters.avoid_otp import AvoidOTP


def avoid_otp(
        sam_path: str,
        out_path: str,
        target_region: str,
        density_thresh: float = 1e-4,
        search_range: t.Tuple[int, int] = (-5e5, 5e5),
        ):
    region = parse_region(target_region)
    blocks = list(read_align_blocks(sam_path))
    blocks.sort(key=lambda b: count_ovlap_with_region(b, region)[0], reverse=True)
    acc = AvoidOTP(region, density_thresh, search_range)
    blocks = acc.filter(blocks)
    counted = []
    for b in blocks:
        c = count_ovlap_with_region(b, region)
        if c[0] > 0:
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
