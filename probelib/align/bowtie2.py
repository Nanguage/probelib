import subprocess as subp


def align_se_sen(fq_path: str,
                 index: str,
                 sam_path: str,
                 threads: int = 10,
                 header: bool = False,
                 ) -> str:
    cmd = ["bowtie2", "-x", index, "-U", fq_path]
    if not header:
        cmd.append("--no-hd")
    cmd += ["-t", "-k", "100", "--very-sensitive-local"]
    cmd += ["-p", str(threads)]
    cmd += ["-S", sam_path]
    subp.check_call(cmd)
    return sam_path
